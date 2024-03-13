// Define area of interest (AOI)
var AOI = ee.Geometry.Polygon( [
  [-65.77070635387585,-17.0787038444975],
  [-65.15097045709443,-17.0787038444975],
  [-65.15097045709443,-16.29066732130244],
  [-65.77070635387585,-16.29066732130244],
  [-65.77070635387585,-17.0787038444975]
]);

// Get two DEM options for terrain flattening 
var dem_srtm = ee.Image('USGS/SRTMGL1_003').clip(AOI);
var dem_cop = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM').filterBounds(AOI).mosaic();

// -------------------Fethcing and preprocessing Sentinel-1 images--------------------//
// Sentinel-1 preprocessing code is from https://github.com/adugnag/gee_s1_ard

var wrapper = require('users/adugnagirma/gee_s1_ard:wrapper');
var helper = require('users/adugnagirma/gee_s1_ard:utilities');

// DEFINE PARAMETERS
var parameter = {//1. Data Selection
    START_DATE: "2020-07-01",
    STOP_DATE: "2020-09-01",
    POLARIZATION:'VVVH',
    ORBIT : 'DESCENDING',
    GEOMETRY: AOI, 
    //2. Additional Border noise correction
    APPLY_ADDITIONAL_BORDER_NOISE_CORRECTION: true,
    //3.Speckle filter
    APPLY_SPECKLE_FILTERING: true,
    SPECKLE_FILTER_FRAMEWORK: 'MULTI',
    SPECKLE_FILTER: 'LEE',
    SPECKLE_FILTER_KERNEL_SIZE: 9,
    SPECKLE_FILTER_NR_OF_IMAGES: 10,
    //4. Radiometric terrain normalization
    APPLY_TERRAIN_FLATTENING: true,
    DEM: dem_cop, //dem_srtm,
    TERRAIN_FLATTENING_MODEL: 'VOLUME',
    TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER: 0,
    //5. Output
    FORMAT : 'DB',
    CLIP_TO_ROI: false,
    SAVE_ASSETS: false
};

// Preprocess the S1 collection
var s1_preprocess = wrapper.s1_preproc(parameter);
// Unprocessed image collection is stored at index 0, preprocessed image collection stored at index 1.
var s1 = s1_preprocess[0];
s1_preprocess = s1_preprocess[1];

// select one image from time series to create function.
var imageS1 = ee.Image(s1_preprocess.median()).select('VV','VH').clip(AOI);

// Visualisation parameters
var s1_band='VH';
var vizParams = {bands: s1_band, min: 0.0032,  max: 0.31};
Map.addLayer(imageS1, vizParams, 'Sentinel-1 VH image',false);
//VV band Visualisation parameters
var s1_band='VV';
var vizParams = {bands: s1_band,min: 0.01, max: 1};
Map.addLayer(imageS1, vizParams, 'Sentinel-1 VV image',false);

// -------------------Fethcing and preprocessing Sentinel-2 images--------------------//
var sentinel2 = ee.ImageCollection("COPERNICUS/S2");
var s2_col = sentinel2.filterBounds(AOI)
              .filterDate('2020-07-01','2020-09-01')
              .filterMetadata('CLOUDY_PIXEL_PERCENTAGE','less_than', 30)
              .map(addVariables);
              
print ('Sentinel-2 Collection',s2_col);

var trueColor_sentinel2 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.2};                  
var s2_image = s2_col.map(maskS2clouds).median();
Map.addLayer(s2_image, trueColor_sentinel2, 'Sentinel-2 True Color');

var imageS2 = addVariables(s2_image).clip(AOI).select(['MNDWI', 'NDWI', 'NDBI', 'NDVI']);

Map.addLayer(imageS2.select('MNDWI'),{min: -0.8, max: 0.8, palette: ['green','white','blue']},'Cloudless Sentinel-2 MNDWI',false);
Map.addLayer(imageS2.select('NDWI'),{min: -1, max: 0.3, palette: ['green','white','purple']},'Cloudless Sentinel-2 NDWI',false);
Map.addLayer(imageS2.select('NDBI'),{min: -1, max: 1, palette: ['cyan', 'red']},'Cloudless Sentinel-2 NDBI',false);
Map.addLayer(imageS2.select('NDVI'), {min: -1, max: 1, palette: ['blue', 'white', 'green']},'Cloudless Sentinel-2 NDVI',false);

//------------------- Normalize all the band pixel values between 0 and 1 -------------//
var image_s1 = normalizeImage(imageS1);
var image_s2 = normalizeImage(imageS2);
var image_s12 = image_s1.addBands(image_s2);
print(image_s12)

Map.centerObject(AOI,13);
Map.addLayer(AOI, {color: 'yellow'}, 'AOI');

//------------------ Collecting training data ------------------//
// This step requires a user add training polygons on the map for 4 classes.
// Assign the 'landcover' property as follows://

// vege: 0
// water: 1
// sedi: 2
// urban_agri: 3

//------------------- Supervised classification ----------------//
//Merge the training data collected from above section into one FeatureCollection
var samples = water.merge(sedi).merge(vege).merge(urban_agri);

var sample_img = ee.Image().byte().paint(samples, "landcover").rename("landcover")

// Then apply stratified sampling to this land cover sample image
var stratifiedsample=sample_img.stratifiedSample({
  numPoints: 1000, 
  classBand: 'landcover',
  scale: 30, 
  projection: 'EPSG:4326',
  region: AOI2,
  geometries:true
});

print('Stratified samples', stratifiedsample); 
print(stratifiedsample.reduceColumns(ee.Reducer.frequencyHistogram(),['landcover']).get('histogram','No of points')); 

// Add a random value field to the sample and use it to approximately split 70%
// of the features into a training set and 20% into a validation set.
var all_sample = stratifiedsample.randomColumn();
var trainSample = all_sample.filter('random <= 0.7');
var validSample = all_sample.filter('random > 0.7');
print(trainSample);
print(validSample);

var bands_s1 = ['VH','VV']; //'logVV/VH'
var bands_s2 = ['NDWI','MNDWI','NDBI','NDVI'];
var bands_s12 = ['VH','VV','NDWI','MNDWI','NDBI','NDVI'];

var training_s1 = image_s1.select(bands_s1).sampleRegions({
collection: trainSample,
properties: ['landcover'],
scale: 30
});
var training_s2 = image_s2.select(bands_s2).sampleRegions({
collection: trainSample,
properties: ['landcover'],
scale: 30
});
var training_s12 = image_s12.select(bands_s12).sampleRegions({
collection: trainSample,
properties: ['landcover'],
scale: 30
});

//Train classifier - e.g. smileCart(), randomForest, svm
var classifier_s1 = ee.Classifier.smileRandomForest(10).train({
features: training_s1,
classProperty: 'landcover',
inputProperties: bands_s1
});
var classifier_s2 = ee.Classifier.smileRandomForest(10).train({
features: training_s2,
classProperty: 'landcover',
inputProperties: bands_s2
});
var classifier_s12 = ee.Classifier.smileRandomForest(10).train({
features: training_s12,
classProperty: 'landcover',
inputProperties: bands_s12
});

//Run the classification
var classified_s1 = image_s1.select(bands_s1).classify(classifier_s1);
var classified_s2 = image_s2.select(bands_s2).classify(classifier_s2);
var classified_s12 = image_s12.select(bands_s12).classify(classifier_s12);

Map.addLayer(classified_s1,
{min: 0, max: 3, palette:  ['green', 'blue', 'orange', 'grey']},'Sentinel-1 classification');
Map.addLayer(classified_s2,
{min: 0, max: 3, palette:  ['green', 'blue', 'orange', 'grey']},'Sentinel-2 classification');
Map.addLayer(classified_s12,
{min: 0, max: 3, palette:  ['green', 'blue', 'orange', 'grey']},'Sentinel-1&2 classification');

// ----------------------------Validate the supervised classification--------------------//
var validation_s1 = classified_s1.sampleRegions({
collection: validSample,
properties: ['landcover'],
scale: 30,
});
print('validation for Sentinel-1:', validation_s1);
var validation_s2 = classified_s2.sampleRegions({
collection: validSample,
properties: ['landcover'],
scale: 30,
});
print('validation for Sentinel-2:', validation_s2);
var validation_s12 = classified_s12.sampleRegions({
collection: validSample,
properties: ['landcover'],
scale: 30,
});
print('validation for Sentinel-1&2:', validation_s12);
//Compare the landcover of your validation data against the classification result
var testAccuracy_s1 = validation_s1.errorMatrix('landcover', 'classification');
print('Sentinel 1 Validation error matrix: ', testAccuracy_s1);
print('Sentinel 1 Validation overall accuracy: ', testAccuracy_s1.accuracy());
print('Sentinel 1 Validation Kappa coeeficient: ', testAccuracy_s1.kappa());

var testAccuracy_s2 = validation_s2.errorMatrix('landcover', 'classification');
print('Sentinel 2 Validation error matrix: ', testAccuracy_s2);
print('Sentinel 2 Validation overall accuracy: ', testAccuracy_s2.accuracy());
print('Sentinel 2 Validation Kappa coeeficient: ', testAccuracy_s2.kappa());

var testAccuracy_s12 = validation_s12.errorMatrix('landcover', 'classification');
print('Sentinel 1&2 Validation error matrix: ', testAccuracy_s12);
print('Sentinel 1&2 Validation overall accuracy: ', testAccuracy_s12.accuracy());
print('Sentinel 1&2 Validation Kappa coeeficient: ', testAccuracy_s12.kappa());

// --------------------- Export classification map and S2 base map --------------------//
var AOI2 = ee.Geometry.Polygon( [
  [-65.76383989879773,-16.87643392893412],
  [-65.37344360162568,-16.87643392893412],
  [-65.37344360162568,-16.49750798029421],
  [-65.76383989879773,-16.49750798029421],
  [-65.76383989879773,-16.87643392893412]
]);

Export.image.toAsset({
image: imageS2,
description: 'S2TrueColorImage',
assetId: 'Bolivia2020_example_S2_trueColor',
crs: projection.crs,
scale: 30,
region: AOI2
});

Export.image.toAsset({
// Export.image.toDrive({
image: classified_s12,
description: 'ClassMap',
assetId: 'Bolivia2020_example_S12_Channel_class',
crs: projection.crs,
scale: 30,
region: AOI2
//   });
});

//----------------------------- FUNCTIONS -----------------------------------//
function addVariables(image){
//method from https://zhuanlan.zhihu.com/p/136929576
var mndwi = image.normalizedDifference(['B3','B11']).rename('MNDWI');
var ndwi = image.normalizedDifference(['B3','B8']).rename('NDWI');
var ndbi = image.normalizedDifference(['B11','B8']).rename('NDBI');
var ndvi = image.normalizedDifference(['B8','B4']).rename('NDVI');
return image.addBands([mndwi,ndwi,ndbi,ndvi]);
}

function maskS2clouds(image) {
var qa = image.select('QA60');

// Bits 10 and 11 are clouds and cirrus, respectively.
var cloudBitMask = 1 << 10;
var cirrusBitMask = 1 << 11;

// Both flags should be set to zero, indicating clear conditions.
var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
.and(qa.bitwiseAnd(cirrusBitMask).eq(0));

return image.updateMask(mask).divide(10000)
.select("B.*")
.copyProperties(image,["system:time_start"]);
}

function normalizeImage(image){
//method from: https://gis.stackexchange.com/questions/313394/normalization-in-google-earth-engine
// calculate the min and max value of an image
var minMax = image.reduceRegion({
reducer: ee.Reducer.minMax(),
geometry: image.geometry(),
scale: 30,
maxPixels: 10e9,
tileScale: 2
}); 
// use unit scale to normalize the pixel values
var unitScale = ee.ImageCollection.fromImages(
image.bandNames().map(function(name){
name = ee.String(name);
var band = image.select(name);
return band.unitScale(ee.Number(minMax.get(name.cat('_min'))), ee.Number(minMax.get(name.cat('_max'))))
        // eventually multiply by 100 to get range 0-100
        //.multiply(100);
})).toBands().rename(image.bandNames());

print(unitScale)
return unitScale
}

function otsu(histogram) {
var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
var size = means.length().get([0]);
var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
var mean = sum.divide(total);

var indices = ee.List.sequence(1, size);

// Compute between sum of squares, where each mean partitions the data.
var bss = indices.map(function(i) {
var aCounts = counts.slice(0, 0, i);
var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
var aMeans = means.slice(0, 0, i);
var aMean = aMeans.multiply(aCounts)
.reduce(ee.Reducer.sum(), [0]).get([0])
.divide(aCount);
var bCount = total.subtract(aCount);
var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
return aCount.multiply(aMean.subtract(mean).pow(2)).add(
bCount.multiply(bMean.subtract(mean).pow(2)));
});

print(ui.Chart.array.values(ee.Array(bss), 0, means));

// Return the mean value corresponding to the maximum BSS.
return means.sort(bss).get([-1]);
}

//Dynamic Otsu thresholding method including Canny edge detection and buffering water pixels
// Source: https://code.earthengine.google.com/e9c699d7b5ef811d0a67c02756473e9d
function computeThresholdUsingOtsu(image, scale, bounds, cannyThreshold, cannySigma, minValue, debug, minEdgeLength, minEdgeGradient, minEdgeValue) {
// clip image edges
var mask = image.mask().gt(0).clip(bounds).focal_min(ee.Number(scale).multiply(3), 'circle', 'meters');

// detect sharp changes
var edge = ee.Algorithms.CannyEdgeDetector(image, cannyThreshold, cannySigma);
edge = edge.multiply(mask);

if(minEdgeLength) {
var connected = edge.mask(edge).lt(cannyThreshold).connectedPixelCount(200, true);

var edgeLong = connected.gte(minEdgeLength);

if(debug) {
print('Edge length: ', ui.Chart.image.histogram(connected, bounds, scale, buckets));

Map.addLayer(edge.mask(edge), {palette:['ff0000']}, 'edges (short)', false);
}

edge = edgeLong;
}

// buffer around NDWI edges
var edgeBuffer = edge.focal_max(ee.Number(scale), 'square', 'meters');

if(minEdgeValue) {
var edgeMin = image.reduceNeighborhood(ee.Reducer.min(), ee.Kernel.circle(ee.Number(scale), 'meters'))

edgeBuffer = edgeBuffer.updateMask(edgeMin.gt(minEdgeValue))

if(debug) {
Map.addLayer(edge.updateMask(edgeBuffer), {palette:['ff0000']}, 'edge min', false);
}
}

if(minEdgeGradient) {
var edgeGradient = image.gradient().abs().reduce(ee.Reducer.max()).updateMask(edgeBuffer.mask())

var edgeGradientTh = ee.Number(edgeGradient.reduceRegion(ee.Reducer.percentile([minEdgeGradient]), bounds, scale).values().get(0))

if(debug) {
print('Edge gradient threshold: ', edgeGradientTh)

Map.addLayer(edgeGradient.mask(edgeGradient), {palette:['ff0000']}, 'edge gradient', false);

print('Edge gradient: ', ui.Chart.image.histogram(edgeGradient, bounds, scale, buckets))
}

edgeBuffer = edgeBuffer.updateMask(edgeGradient.gt(edgeGradientTh))
}

edge = edge.updateMask(edgeBuffer)
var edgeBuffer = edge.focal_max(ee.Number(scale).multiply(1), 'square', 'meters');
var imageEdge = image.mask(edgeBuffer);

if(debug) {
Map.addLayer(imageEdge, {palette:['222200', 'ffff00']}, 'image edge buffer', false)
}

// compute threshold using Otsu thresholding
var buckets = 100;
var hist = ee.Dictionary(ee.Dictionary(imageEdge.reduceRegion(ee.Reducer.histogram(buckets), bounds, scale)).values().get(0));

var threshold = ee.Algorithms.If(hist.contains('bucketMeans'), otsu(hist), minValue);
threshold = ee.Number(threshold)


if(debug) {
// experimental
// var jrc = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('occurrence')
// var jrcTh = ee.Number(ee.Dictionary(jrc.updateMask(edge).reduceRegion(ee.Reducer.mode(), bounds, scale)).values().get(0))
// var water = jrc.gt(jrcTh)
// Map.addLayer(jrc, {palette: ['000000', 'ffff00']}, 'JRC')
// print('JRC occurrence (edge)', ui.Chart.image.histogram(jrc.updateMask(edge), bounds, scale, buckets))

Map.addLayer(edge.mask(edge), {palette:['ff0000']}, 'edges', true);

print('Threshold: ', threshold);

print('Image values:', ui.Chart.image.histogram(image, bounds, scale, buckets));
print('Image values (edge): ', ui.Chart.image.histogram(imageEdge, bounds, scale, buckets));
Map.addLayer(mask.mask(mask), {palette:['000000']}, 'image mask', false);
}

return minValue !== 'undefined' ? threshold.max(minValue) : threshold;
}
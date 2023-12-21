var AOI2 = ee.Geometry.Polygon( [
    [-65.76383989879773,-16.87643392893412],
    [-65.37344360162568,-16.87643392893412],
    [-65.37344360162568,-16.49750798029421],
    [-65.76383989879773,-16.49750798029421],
    [-65.76383989879773,-16.87643392893412]
  ]);

//--------------- Load the classification map and S2 base map you generated using GEE_supervised_classification_Sentinel_1_2.js ---------//

// Change the content in the " " to the asset name that you customized for your own
var classified = ee.Image("users/qiuyangschen/Bolivia2020_example_S12_Channel_class");
var S2_trueColor = ee.Image("users/qiuyangschen/Bolivia2020_example_S2_trueColor");

Map.addLayer(S2_trueColor,{bands: ['B4', 'B3', 'B2'], min: 0, max: 0.2},'Sentinel2 B432 Bolivia2020');
Map.addLayer(classified,
{min: 0, max: 3, palette:  ['green', 'blue', 'orange','grey']},'Sentinel-1&2 classification',false); 

//----------------------Remove pepper-and-salt errors: Clustering ----------------//
var seeds = ee.Algorithms.Image.Segmentation.seedGrid(5);

var snic = ee.Algorithms.Image.Segmentation.SNIC({
  image: S2_trueColor.select('B.*'), 
  compactness: 0,
  connectivity: 8,
  neighborhoodSize: 50,
  size: 2,
  seeds: seeds
})
var clusters = snic.select('clusters')
// Assign class to each cluster based on 'majority' voting (using ee.Reducer.mode()
var smoothed = classified.addBands(clusters);

var clusterMajority = smoothed.reduceConnectedComponents({
  reducer: ee.Reducer.mode(),
  labelBand: 'clusters'
});
Map.addLayer(clusterMajority, {min: 0, max: 3, palette: ['blue', 'orange','green','grey']},'Classmap after clustering');

//--------------------------Export images------------------------------
var projection = S2_trueColor.select('B2').projection().getInfo();
Export.image.toDrive({
  image:  S2_trueColor.visualize({bands: ['B4', 'B3', 'B2'], min: 0, max: 0.2}),
  description: 'Bolivia2020_example_S2_TrueColor',
  crs: projection.crs,
  scale: 30,
  region: AOI2
  });
  
Export.image.toDrive({
  image: clusterMajority,
  description: 'Bolivia2020_example_Clustered_ClassMap',
  crs: projection.crs,
  scale: 30,
  region: exportAOI
  });

# -*- coding: utf-8 -*-

from sys import exit
from math import hypot, sqrt，pi
from numpy import linspace, zeros
from geopandas import read_file
from numpy.random import randint
from rasterio import open as rio_open
from rasterio.plot import show as rio_show
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import subplots, savefig, get_cmap
from rtree import index
import time
t1 = time.time()


def WeightedRedistribution(w, s, pointData, weightingSurface, polygons):
    '''
    * Main function to achieve the redistribution
    * w is the number of user-designed random points in each polygon.
    ** w - the higher,the more accuary (20)
    * s is the user-designed ambiguity of output map
    ** s - the higher, the more ambiguity (0.1)
    '''
    # get polygons from the gm-district dataset
    for id, poly in polygons.iterrows():

        # For each district, get the tweets that are within it
        # get all the tweets who within admins by rtree
        pointinpolyindex = list(idx.intersection(poly.geometry.bounds))

        pointinpoly = tweets.iloc[pointinpolyindex]

        points = pointinpoly[pointinpoly.within(poly.geometry)]

        # For each tweet, use loop to create new random locations
        for id, p in points.iterrows():

             # get the new location with the highest value on the population surface
            weightedpoint = relocate(p, poly, w, weightingSurface)
            
            # distribution radius
            r = radius(poly.geometry.area, s)
            
            # calculate cell value in distribution
            # use list comprehension as a shorthand way of creating a nested loop
            distribution = [(i, j) for i in range(weightedpoint[0] - r, weightedpoint[0] + r) for j in range(weightedpoint[1] - r, weightedpoint[1] + r)]
            
            for i in distribution:
                
                if computeDistance(weightedpoint[0], weightedpoint[1], i[0], i[1]) <= r :
                    
                    # add your weighted distribution of values to the output raster dataset
                    output[i] += cellValue(weightedpoint[0], weightedpoint[1], i[0], i[1], r)
    return output


def relocate(point, poly, iterations, weightingSurface):
    '''
    * Relocates a point using the weighting surface.
    * the more random points the more accuracy
    * renturns the hightest value of weightedsurface in image points
    '''

    print(point.geometry)
    
    # create random points 
    randomX = randint(poly.geometry.bounds[0], poly.geometry.bounds[2], size=iterations)
    
    randomY = randint(poly.geometry.bounds[1], poly.geometry.bounds[3], size=iterations)

    # holds max value
    maxVal = 0
    weightedpoint = 0
     
    # create the required number of offset points
    # offset the point and push into an array list
    for x, y in zip(randomX, randomY):
        # convert coord to rowcol
        weightedimage = coord2Image(pop, x, y)
        
        # extract the value from pop
        value = pop_data[weightedimage]
        
        # find the hightest value
        if value > maxVal:
            
            maxVal = value
            
            weightedpoint = weightedimage
                
    return weightedpoint


def computeDistance(x1, y1, x2, y2):
    '''
    * calculate the the Euclidean norm
    * returns the distance
    '''
    d = hypot(x2 - x1, y2 - y1)

    return d


def cellValue(x1, y1, x2, y2, r):
    '''
    * calculate cell value in distribution around seed
    * returns the value
    '''
    v = 1 - sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2) / r

    return v


def radius(area, s):
    '''
    * calculate the distrubution radius based on the area of polygon and user-set Ambiguity
    * returns the int r in coords 
    '''
    r = sqrt(area * s / pi)
    
    return int(r/resolution)


def coord2Image(d, x, y):
    '''
    * convert between coords and array position
    * returns row,col (y,x) as expected by rasterio
    '''
    r, c = d.index(x, y)
    return int(r), int(c)


# It's as if the interpreter inserts this at the top
# of your module when run as the main program.
if __name__ == "__main__":
    
    '''
    * read three input datasets
    ''' 
    # administrative areas - districts
    administrativeAreas = read_file("./wr/gm-districts.shp")

    # level3
    lsoa = read_file("./wr/gm-lsoa.shp")

    # tweets points to be redistributed
    tweets = read_file("./wr/level3-tweets-subset.shp")

    '''
    * generate weight redistribution
    '''
    # initialise an rtree index
    idx = index.Index()

    # loop through each row and construct spatial index on the larger dataset (population)
    for id, tw in tweets.iterrows():
        idx.insert(id, tw.geometry.bounds)

    # pop as the weighting surface and process the weighted rdistribution
    with rio_open("./wr/GM_POP.tif") as pop:

        # print(pop.profile)
        # get value from weighting surface at p[i]
        # get pop value
        pop_num = pop.read()[0]

        # get resolution
        resolution = pop.res[0]
        # print(resolution)

        # convert meter to pixel value
        pop_data = pop_num/resolution

        # create a output dataset of 0's that matches the dimensions of pop
        output = zeros((pop.height, pop.width))
        
        # generate the weight redistrubution and catch exception
        try:
            WR_map = WeightedRedistribution(20, 0.1, tweets, pop, administrativeAreas)
            
        except IndexError:
            print("Failed, Please increase w or decrease s to fit the bound")
            exit()
    
    '''
    # Plot the output surface
    '''        
    # plot setting
    fig, ax = subplots(1, figsize=(12, 12))
    ax.set_title('The Population-weighted Redistribution',fontsize=25)
    
    # remove axes
    ax.axis('off')
    
    # set map palette
    cMap = get_cmap('RdPu')
    
    # add arrow to map
    x, y, arrow_length = 0.05, 0.9, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
            arrowprops=dict(facecolor='black', width=5, headwidth=15),
            ha='center', va='center', fontsize=20,
            xycoords=ax.transAxes)
    
    # add raster colorbar to map
    img = ax.imshow(WR_map, cmap=cMap)
    position = fig.add_axes([0.15,0.23,0.02,0.2])
    fig.colorbar(img, cax = position, label = 'Tweet Activity')
    
    # add reference layer to map
    '''
    administrativeAreas.plot(
        ax = ax,
        facecolor="none",
        edgecolor='black',
        linewidth = 0.7
        )
    '''
    lsoa.plot(
        ax=ax,
        edgecolor="grey",
        facecolor="none",
        linewidth=0.7,
    )
    
    tweets.plot(
        ax=ax,
        color="Blue"
    )
    
    # add legend to map
    ax.legend(labels=['Original Tweets'], fontsize=20, loc='lower left')
    
    # plot raster layer to map
    rio_show(WR_map, ax =ax, transform = pop.transform, cmap = LinearSegmentedColormap.from_list('RePu', cMap(linspace(0, 1, 100))))

    # save the result
    savefig('./GM_WR_RW.png')
    print(f'{time.time() - t1}')
    print("done!")

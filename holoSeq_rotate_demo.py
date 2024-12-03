# see https://github.com/fubar2/holoSeq
# illustrates some of the basic ideas in converting
# a set of features that have been mapped to a genome into a
# linear or in this case 2D display.
# 
# The diagonal is rotated to become the horizontal axis
# and the tap coordinates are back calculated to the original unrotated ones 
# so the genomic position on the original diagonal 2D plot are shown
# 
# run with
# panel serve holoSeq_rotate_demo.py --show
# should pop up 10000 points.
# try 10 million - still works smoothly
# Needs depenencies
# pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh
#

from bisect import bisect_left
from collections import OrderedDict
import gzip
import itertools
import math
import numpy as np

import holoviews as hv
import pandas as pd
import panel as pn

from rotater import rotater

from holoviews.operation.datashader import (
    rasterize,
    dynspread,
)

hv.extension("bokeh")
pn.extension()


SQRT2 = math.sqrt(2)

def showTap(x, y, rotated, rot):
    # populate a string widget with tap coordinates translated into contig and offset
    if np.isnan(x) or np.isnan(y):
        s = "Click on plot for coordinates as contig:offset"
    else:
        chrx = "Out of range"
        offsx = 0
        chry = "Out of range"
        offsy = 0
        xur =  yur = None
        if rotated:
            xur, yur = rot.unrotatecoords(xr=x, yr=y)
            i = bisect_left(hstarts, xur)
            if i > 0 and i <= len(hnames):
                chrx = hnames[i - 1]
                offsx = xur - hstarts[i - 1]
            i = bisect_left(hstarts, yur)
            if i > 0 and i <= len(hnames):
                chry = hnames[i - 1]
                offsy = yur - hstarts[i - 1]
            s = "Rotated X axis genome %s:%d Rotated Y axis genome %s:%d x %d y %d xur %d yur %d" % (
            chrx,
            offsx,
            chry,
            offsy,
            x,
            y,
            xur,
            yur)
        else:
            i = bisect_left(hstarts, x)
            if i > 0 and i <= len(hnames):
                chrx = hnames[i - 1]
                offsx = x - hstarts[i - 1]
            i = bisect_left(hstarts, y)
            if i > 0 and i <= len(hnames):
                chry = hnames[i - 1]
                offsy = y - hstarts[i - 1]
            s = "X axis %s:%d Y axis %s:%d. x %d y %d" % (chrx, offsx, chry, offsy, x, y)
    str_pane = pn.pane.Str(
        s,
        styles={"font-size": "10pt", "color": "darkblue", "text-align": "center"},
        width=width,
    )
    return str_pane

SQRT2 = math.sqrt(2)
xwidth = 3000000
xmax = 10000
width = 1000
height = 400
rng = np.random.default_rng(1)  # all plots will be identical !
hlen = [xwidth / 2, xwidth / 4, xwidth / 8, xwidth / 8]
hstarts = list(itertools.accumulate(hlen))
hstarts.insert(0, 0)
hnames = ["chr%d" % i for i in range(1, 5)]
title = "%d random points" % xmax
ticks = [(hstarts[i], hnames[i]) for i in range(4)]
xcoords = np.array([rng.uniform(0, xwidth) for i in range(xmax)])
ycoords = np.array([xcoords[i] - abs(rng.uniform(0, xcoords[i])) for i in range(xmax)])
points = hv.Points((xcoords, ycoords))
stream = hv.streams.Tap(source=points, x=np.nan, y=np.nan)
rot = rotater(xwidth, xwidth)
showloc = pn.bind(showTap, x=stream.param.x, y=stream.param.y, rotated=False, rot=rot)
xyr = rot.rotatecoords(
    xcoords,
    ycoords,
    adjust=True
)
rotpoints = hv.Points(xyr, kdims=["x", "y"])
streamr = hv.streams.Tap(source=rotpoints, x=np.nan, y=np.nan)
showlocr = pn.bind(showTap, x=streamr.param.x, y=streamr.param.y, rotated=True, rot=rot)
rot = pn.pane.HoloViews(
    dynspread(
        rasterize(rotpoints)
        .relabel("%s" % title)
        .opts(
            # ylim=(0,xwidth),
            framewise=True,
            autorange=None,
            cmap="inferno",
            cnorm="log",
            colorbar=True,
            width=width,
            height=height,
            xticks=ticks,
            yticks=[(0,0)],
            xrotation=45,
            fontsize={"xticks": 7, "yticks": 7},
            scalebar=True,
            scalebar_range="x",
            scalebar_location="center_left",
            scalebar_unit=("bp"),
            shared_axes=False,
            tools=["xwheel_zoom",
                "box_zoom",
                "tap",
                "pan",
                "reset",
            ],
            default_tools=[],
            active_tools=["xwheel_zoom", "tap", "pan"],
        )
    )
)

unrot = pn.pane.HoloViews(
    dynspread(
        rasterize(points)
        .relabel("%s" % title)
        .opts(
            cmap="inferno",
            cnorm="log",
            colorbar=True,
            width=width,
            height=height,
            xticks=ticks,
            yticks=ticks,
            xrotation=45,
            fontsize={"xticks": 7, "yticks": 7},
            tools=["tap"],
            scalebar=True,
            scalebar_range="x",
            scalebar_location="top_left",
            scalebar_unit=("bp"),
            shared_axes=False,
        )
    )
)
pnc = pn.Column(showlocr, rot, showloc, unrot)


pnc.servable(
    title=title,
)
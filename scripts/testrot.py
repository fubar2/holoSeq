import math
import pandas as pd
import numpy as np
import random

SQRT2 = math.sqrt(2)


class rotater:
    """
    moved into a class so the constants are only calculated once
    if single pairs are being rotated
    x 0 y 0 xr -2121320 yr 878679
    x 3000000 y 0 xr 0 yr -1242640
    x 1500000 y 1500000 xr 0 yr 878679
    x 3000000 y 3000000 xr 2121320 yr 878679
    """

    def __init__(self, xwidth, ywidth):
        self.rotate = True
        self.onepointRot = True
        self.radians = 0.7853981633974483
        self.origin = (0, ywidth)
        self.cos_rad = math.cos(self.radians)
        self.sin_rad = math.sin(self.radians)
        (self.xmin, self.ymax) = self.rotatecoords(0, 0, adjust=False)
        self.ymin = self.rotatecoords(xwidth, 0, adjust=False)[1]
        self.xmax = self.rotatecoords(xwidth, ywidth, adjust=False)[0]
        self.xnew = SQRT2 * xwidth
        self.ynew = ywidth / SQRT2
        self.xwidth = xwidth
        self.ywidth = ywidth
        self.xscalefact = self.xwidth / self.xnew
        self.yscalefact = self.ywidth / self.ynew
        print(
            "xmin %d xnew %d ymin %d ynew %d, cos %f, sin %f"
            % (self.xmin, self.xnew, self.ymin, self.ynew, self.cos_rad, self.sin_rad)
        )

    def rotatecoords(
        self,
        xin=0,
        yin=0,
        adjust=True,
    ):
        """
        This version optimised to operate on pair at a time so precomputed constants
        # make a rotated heatmap where the diagonal becomes the horizontal x axis
        # https://gist.github.com/LyleScott/d17e9d314fbe6fc29767d8c5c029c362
        # this inflates the xaxis by sqrt(2) but can rescale and it seems to work (TM)
        # xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
        # pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
        # y height is complicated - h^2 + (1/2*sqrt2*xwdith)^2 = y^2
        """
        self.offset_x, self.offset_y = self.origin
        adjusted_x = xin - self.offset_x
        adjusted_y = yin - self.offset_y
        qx = self.offset_x + self.cos_rad * adjusted_x + self.sin_rad * adjusted_y
        qy = self.offset_y + -self.sin_rad * adjusted_x + self.cos_rad * adjusted_y
        if adjust:
            xrs = (qx - self.xmin) * self.xscalefact
            yrs = (qy - self.ymin) * self.yscalefact
        else:
            xrs = qx
            yrs = qy
        if self.onepointRot:
            return (xrs, yrs)
        else:
            xyr = pd.DataFrame(np.vstack([xrs, yrs]).T, columns=["x", "y"])
            return xyr

    def unrotatecoords(self, xr, yr):
        """
        rotated coords back calculate
        """
        qx = xr / self.xscalefact + self.xmin
        qy = yr / self.yscalefact + self.ymin
        adjx = (self.offset_y - qy + qx - self.offset_x) / 2
        adjy = qx - self.offset_x - adjx
        x = adjx / self.cos_rad + self.offset_x
        y = adjy / self.cos_rad + self.offset_y
        return (x, y)


xwidth = 3000000
r = rotater(xwidth, xwidth)
r.onepointRot = True
tests = [(0, 0), (xwidth, 0), (xwidth / 2, xwidth / 2), (xwidth, xwidth)]
if False:
    for x, y in tests:
        xr, yr = r.rotatecoords(x, y)
        xrs = (xr - r.xmin) * r.xscalefact
        yrs = (yr - r.ymin) * r.yscalefact
        print("x %d y %d xr %d yr %d xrs %d yrs %d" % (x, y, xr, yr, xrs, yrs))
else:
    for i in range(50):
        y = random.randint(0, xwidth)
        x = max(0, y - random.randint(0, 100))
        xr, yr = r.rotatecoords(x, y)
        xur, yur = r.unrotatecoords(xr, yr)
        print("x %d y %d xr %d yr %d xur %d yur %d" % (x, y, xr, yr, xur, yur))

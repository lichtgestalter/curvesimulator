"""
There are n circles with known position and radius. Some circles overlap.
Calculate the visible relative area of each circle.
Complete the function visible_percentages().
Consider the funtion intersection() to be correctly implemented.
"""

from shapely.geometry import Point  # https://shapely.readthedocs.io/en/stable/manual.html
import math

class TCircle():

    def __init__(self, x, y, r):
        self.x = x  # x-coordinate
        self.y = y  # y-coordinate
        self.r = r  # radius


# Example circles
eclipsee = TCircle(10, 20, math.sqrt(1/math.pi))
eclipser1 = TCircle(10, 21.79, 0.567)
eclipser2 = TCircle(10, 20.9, 0.567)

circle_list = [eclipsee, eclipser1, eclipser2]

shapely_circles = [Point(c.x, c.y).buffer(c.r, resolution=4096) for c in circle_list]  # Create shapely circles
intersection01 = shapely_circles[0].intersection(shapely_circles[1])
intersection02 = shapely_circles[0].intersection(shapely_circles[2])
intersection012 = shapely_circles[0].intersection(shapely_circles[1]).intersection(shapely_circles[2])
eclipsed_rel_0 = (intersection01.area + intersection02.area - intersection012.area) / shapely_circles[0].area


print(f"{shapely_circles[0].area=} {intersection01.area=} {intersection02.area=} {intersection012.area=} {eclipsed_rel_0=}")

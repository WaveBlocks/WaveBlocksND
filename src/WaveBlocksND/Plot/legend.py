"""The WaveBlocks Project

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011 R. Bourquin
@license: Modified BSD License
"""

import pylab


def legend(*args, **kwargs):
    """Overwrites the pylab legend function.

    It adds another location identfier 'outer right'
    which locates the legend on the right side of the plot

    The args and kwargs are forwarded to the pylab legend function
    """
    if kwargs.has_key('loc'):
        loc = kwargs['loc']

        if loc == "outer right":
            kwargs.pop('loc')

            axes = pylab.gca()
            bbox = axes.get_position().get_points()

            x = bbox[0][0]
            y = bbox[0][1]
            w = bbox[1][0]
            h = bbox[1][1]
            
            axes.set_position([x,y,0.5,0.5])

            #legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  
            leg = pylab.legend(loc=(0,0), *args, **kwargs)
            frame = leg.get_frame()
            lw = frame.get_width()
            lh = frame.get_height()
            
            scale = lw / (1.0*w)
            
            leg._loc = (x+w,y+h-lh)

            axes.set_position([x,y,0.7*((w)/scale),0.9*h])

            pylab.draw_if_interactive()
            return leg

    return pylab.legend(*args, **kwargs)

from mvq import *

def t(res, exp):

    if res != exp:
        raise "FAILED"

    print "PASS"

print "bsphere"

bs = bsphere(0, 0, 0, 1)

t(bs.includes(0,0,0), True)
t(bs.includes(vec3d(0,0,0)), True)
t(bs.includes(1,1,1), False)
t(bs.includes(vec3d(1,1,1)), False)
t(bs.includes(bsphere(0,0,0,0.5)), True)
t(bs.includes(bsphere(0,0,0,2)), False)

t(bs.intersects(bsphere(0,0,0,0.5)), True)
t(bs.intersects(bsphere(0,0,0,2)), True)
t(bs.intersects(bsphere(1,1,1,1)), True)

t(bs.intersects(bs), True)
t(bs.includes(bs), True)

print "bbox"

bb = bbox(-1, -1, -1, 1, 1, 1)

t(bb.includes(0,0,0), True)
t(bb.includes(vec3d(0,0,0)), True)
t(bb.includes(1,1,1), True)
t(bb.includes(vec3d(1,1,1)), True)
t(bb.includes(bsphere(0,0,0,0.5)), True)
t(bb.includes(bsphere(0,0,0,2)), False)

t(bb.intersects(bb), True)
t(bb.includes(bb), True)

#t(bb.intersects(bsphere(2, 2, 0, 2)), False)
#t(bb.intersects(bsphere(2, 2, 0, 3)), True)



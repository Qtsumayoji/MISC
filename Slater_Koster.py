from sympy import *

#var("a:z")

ppsigma = Symbol('ppsigma')
pppi = Symbol('pppi')
pdsigma = Symbol('pdsigma')
pdpi = Symbol('pdpi')


theta = Symbol('theta')
phi = Symbol('phi')

ex = Symbol('ex')
ex = sin(theta)*cos(phi)

ey = Symbol('ey')
ey = sin(theta)*sin(phi)

ez = Symbol('ez')
ez = cos(theta)

tx_x = -( ex**2*ppsigma + (1 - ex**2)*pppi )
tx_y = -( ex*ey*ppsigma - ex*ey*pppi )
tx_z = -( ex*ez*ppsigma - ex*ez*pppi )

tx_xy = -( sqrt(3)*ex**2*ey*pdsigma + ey*(1 - 2*ex**2)*pdpi)
tx_yz = -( sqrt(3)*ex*ey*ez*pdsigma - 2*ex*ey*ez*pdpi)
tx_zx = -( sqrt(3)*ex**2*ez*pdsigma + ez*(1 - 2*ex**2)*pdpi)

tx_x2y2 = -( sqrt(3)/2*ex*(ex**2 - ey**2)*pdsigma + ex*(1 - ex**2 + ey**2)*pdpi )
ty_x2y2 = -( sqrt(3)/2*ey*(ex**2 - ey**2)*pdsigma - ey*(1 + ex**2 - ey**2)*pdpi )
tz_x2y2 = -( sqrt(3)/2*ez*(ex**2 - ey**2)*pdsigma - ez*(ex**2 - ey**2)*pdpi )

tx_3z2r2 = -( ex*(ez**2 - (ex**2 + ey**2)/2)*pdsigma - sqrt(3)*ex*ez**2 )
ty_3z2r2 = -( ey*(ez**2 - (ex**2 + ey**2)/2)*pdsigma - sqrt(3)*ey*ez**2 )
tz_3z2r2 = -( ex*(ez**2 - (ex**2 + ey**2)/2)*pdsigma + sqrt(3)*ez*(ex**2 + ey**2) )

def main():
    ppsigma = 0.5
    pppi = -0.3*ppsigma

    th = pi/2
    ph = pi/4
    print("[px-px]")
    print(tx_x.subs([(theta, th), (phi, ph)]),"\n")
    print("[px-py]")
    print(tx_y.subs([(theta, th), (phi, ph)]),"\n")
    print("[px-pz]")
    print(tx_z.subs([(theta, th), (phi, ph)]),"\n")

    th = pi/4
    ph = pi/2
    print("[pz-pz]")
    print(tx_x.subs([(theta, th), (phi, ph)]),"\n")

    th = pi/2
    ph = 0
    print("[px-dx2y2]")
    print(tx_x2y2.subs([(theta, th), (phi, ph)]),"\n")
    print("[py-dx2y2]")
    print(ty_x2y2.subs([(theta, th), (phi, ph)]),"\n")
    print("[pz-dx2y2]")
    print(tz_x2y2.subs([(theta, th), (phi, ph)]),"\n")

    th = pi/2
    ph = 0
    print("[px-3z2r2]")
    print(tx_3z2r2.subs([(theta, th), (phi, ph)]),"\n")
    print("[py-3z2r2]")
    print(ty_3z2r2.subs([(theta, th), (phi, ph)]),"\n")
    print("[pz-3z2r2]")
    print(tz_3z2r2.subs([(theta, th), (phi, ph)]),"\n")


if __name__ == "__main__":
    main()
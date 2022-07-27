

from astropy import units as u
from astropy.coordinates import SkyCoord

radec = ('+0:14:02.86', '-30:22:18.7')
radec = ('+0:13:59.76', '−30:19:29.1')

c = SkyCoord(*radec, unit=(u.hourangle, u.deg))

print(f'{c.ra.deg:.7f} {c.dec.deg:.7f}')

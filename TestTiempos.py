import datetime
from astropy.time import Time

print(datetime.datetime.now())
#print(datetime.datetime.now().jd)
ahora = datetime.datetime.now()
tiempo = Time(ahora, format='datetime', scale='utc')
print(tiempo)
print(tiempo.jd)

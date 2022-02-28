import tools as t


time_jd_utc = t.calendar2julian(2022, 2, 9, 5, 59, 0)
print(time_jd_utc)

time_jd_ut1 = time_jd_utc+(-0.109571/(24*3600))
print(time_jd_ut1)

time_UT1 = (time_jd_ut1 - 2451545)/36525

print(time_UT1)

theta_gmst = 67310.54841 + (876600*3600 + 8640184.812866)*time_UT1 + 0.093104 * time_UT1**2 - 6.2e-6 * time_UT1**3
print(theta_gmst)
theta_gmst = 1/240*(theta_gmst - (theta_gmst//86400 * 86400.0))
print(theta_gmst)

print(t.calendar2julian(2024, 5, 31, 15, 19, 42))
print(t.calendar2julian(2024, 5, 31, 15, 19, 42)-2400000.5)

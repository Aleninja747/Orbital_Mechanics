delta_UT1 = -0.110914
delta_UTC = 37
delta_TT = 32.184

# Dictionary to convert month numbers to names
month_dic = {
    1: 'January',
    2: 'February',
    3: 'March',
    4: 'April',
    5: 'May',
    6: 'June',
    7: 'July',
    8: 'August',
    9: 'September',
    10: 'October',
    11: 'November',
    12: 'December'
}

month_days = {
    1: 31,
    2: 28,
    3: 31,
    4: 30,
    5: 31,
    6: 30,
    7: 31,
    8: 31,
    9: 30,
    10: 31,
    11: 30,
    12: 31
}


def leap_year(year):
    if year % 4 != 0:
        return False
    else:
        if year % 100 == 0:
            return False
        else:
            if year % 400 == 0:
                return True
    return False


def days_in_month(year, month):
    temp_month = month_days
    if leap_year(year):
        temp_month[2] += 1
    return temp_month[month]


class Time:
    def __init__(self, time_sys, year=None, month=None, day=None, hour=None, minute=None, second=None):
        # Initialize all parameters of a time object
        self.__time_sys__ = time_sys
        self.__year__ = year
        self.__month__ = month
        self.__day__ = day
        self.__hour__ = hour
        self.__minute__ = minute
        self.__second__ = second

        # If the year is not inputted we assume the date is in Julian representation
        if year is None:
            self.__representation__ = 'Julian'
            # internally we store the date in modified format to enhance precision
            if day > 2400000.5:
                self.__day__ -= 2400000.5
            self.__hour_representation__ = None

        # If the year is inputted and month is inputted we assume the date is in gregorian representation
        if year is not None and month is not None:
            self.__representation__ = 'Gregorian'

        # If the year is inputted and no month is inputted we assume the date is in date of year representation
        if year is not None and month is None:
            self.__representation__ = 'DOY'

        # If the year is inputted and no hour is inputted we assume the date is in fraction of day representation
        if year is not None and hour is None:
            self.__hour_representation__ = 'FOD'

        # From this point on, the 'day' representation is assumed from the previous checks
        # If the hour is inputted and no minutes are inputted, we assume the time was inputted as a fraction of hour
        if year is not None and hour is not None and minute is None:
            self.__hour_representation__ = 'FOH'
        # If a value is inputted for
        if year is not None and minute is not None:
            self.__hour_representation__ = 'Normal'

        self.time_in_bounds()

    # Overriding string method to allow proper formatting
    def __str__(self):
        if self.__representation__ == 'Julian':
            return str(self.__day__)
        elif self.__representation__ == 'Gregorian' and self.__hour_representation__ == 'Normal':
            return str(self.__year__) + ',' + month_dic[self.__month__] + ',' + str(self.__day__) + ',' + \
                   str(self.__hour__) + ':' + str("%.0f" % self.__minute__) + ':' + str(self.__second__)
        elif self.__representation__ == 'Gregorian' and self.__hour_representation__ == 'FOH':
            return str(self.__year__) + ',' + month_dic[self.__month__] + ',' + str(self.__day__) + ',' + \
                   str(self.__hour__)
        elif self.__hour_representation__ == 'FOD' and self.__representation__ == 'Gregorian':
            return str(self.__year__) + ',' + month_dic[self.__month__] + ',' + str(self.__day__)
        elif self.__representation__ == 'DOY' and self.__hour_representation__ == 'Normal':
            return str(self.__year__) + ',' + str(self.__day__) + ',' + str(self.__hour__) + ':' + \
                   str("%.0f" % self.__minute__) + ':' + str(self.__second__)
        elif self.__representation__ == 'DOY' and self.__hour_representation__ == 'FOD':
            return str(self.__year__) + ',' + str(self.__day__)
        elif self.__representation__ == 'DOY' and self.__hour_representation__ == 'FOH':
            return str(self.__year__) + ',' + str(self.__day__) + ',' + str(self.__hour__)

    # Function that Returns julian day of current time
    def julian(self, modified=True):
        day_term = 0
        hour_term = 0
        if self.__representation__ == 'Julian':
            return self.__day__ + modified * 2400000.5
        elif self.__representation__ == 'DOY':
            day_term = 367 * self.__year__ - int(7 * self.__year__ / 4) + self.__day__
        elif self.__hour_representation__ == 'FOD':
            day_term = 367 * self.__year__ - int((7 * (self.__year__ + int((self.__month__ + 9) / 12))) / 4) + int(
                275 * self.__month__ / 9) + self.__day__ + 1721013.5
        elif self.__representation__ == 'Gregorian':
            day_term = 367 * self.__year__ - int((7 * (self.__year__ + int((self.__month__ + 9) / 12))) / 4) + int(
                275 * self.__month__ / 9) + self.__day__ + 1721013.5
        if self.__hour_representation__ == 'FOH':
            hour_term = self.__hour__ - modified * 2400000.5
        elif self.__hour_representation__ == 'Normal':
            hour_term = + 1 / 24 * (self.__hour__ + (1 / 60) * (self.__minute__ + self.__second__ / 60))
        return day_term + hour_term - modified * 2400000.5

    # Void function that changes the current time system to Terrestrial Time
    def to_tt(self):
        time_diff = 0

        # Set difference between current time system and TT
        if self.__time_sys__ == 'TAI':
            time_diff = delta_TT
        elif self.__time_sys__ == 'UTC':
            time_diff = delta_UTC + delta_TT
        elif self.__time_sys__ == 'UT1':
            time_diff = delta_TT + delta_UTC - delta_UT1

        # Add the time difference to the current time representation
        if self.__representation__ == 'Julian' or self.__hour_representation__ == 'FOD':
            self.__day__ += time_diff / 86400
        elif self.__hour_representation__ == 'FOH':
            self.__hour__ += time_diff / 3600
        else:
            self.__second__ += time_diff
        self.__time_sys__ = 'TT'
        self.time_in_bounds()

    # Void function that changes the current time system to atomic time
    def to_tai(self):
        time_diff = 0

        # Set difference between current time system and TT
        if self.__time_sys__ == 'TT':
            time_diff = -delta_TT
        elif self.__time_sys__ == 'UTC':
            time_diff = delta_UTC
        elif self.__time_sys__ == 'UT1':
            time_diff = delta_UTC + delta_UT1

        # Add the time difference to the current time representation
        if self.__representation__ == 'Julian' or self.__hour_representation__ == 'FOD':
            self.__day__ += time_diff / 86400

        elif self.__hour_representation__ == 'FOH':
            self.__hour__ += time_diff / 3600
        else:
            self.__second__ += time_diff
        self.__time_sys__ = 'TAI'
        self.time_in_bounds()

    # Void function that changes the current time system to UTC
    def to_utc(self):
        time_diff = 0

        # Set difference between current time system and TT
        if self.__time_sys__ == 'TT':
            time_diff = -(delta_TT + delta_UTC)
        elif self.__time_sys__ == 'TAI':
            time_diff = -delta_UTC
        elif self.__time_sys__ == 'UT1':
            time_diff = -delta_UT1

        # Add the time difference to the current time representation
        if self.__representation__ == 'Julian' or self.__hour_representation__ == 'FOD':
            self.__day__ += time_diff / 86400
        elif self.__hour_representation__ == 'FOH':
            self.__hour__ += time_diff / 3600
        else:
            self.__second__ += time_diff
        self.__time_sys__ = 'UTC'
        self.time_in_bounds()

    # Void function that changes the current time system to UT1
    def to_ut1(self):
        time_diff = 0

        # Set difference between current time system and TT
        if self.__time_sys__ == 'TT':
            time_diff = -(delta_TT + delta_UTC) + delta_UT1
        elif self.__time_sys__ == 'TAI':
            time_diff = -delta_UTC + delta_UT1
        elif self.__time_sys__ == 'UTC':
            time_diff = delta_UT1

        # Add the time difference to the current time representation
        if self.__representation__ == 'Julian' or self.__hour_representation__ == 'FOD':
            self.__day__ += time_diff / 86400
        elif self.__hour_representation__ == 'FOH':
            self.__hour__ += time_diff / 3600
        else:
            self.__second__ += time_diff
        self.__time_sys__ = 'UT1'
        self.time_in_bounds()

    def julian_cent(self):
        return (self.julian() - 51544.5) / 36525

    def time_in_bounds(self):
        if self.__representation__ != 'Julian':
            # Start Checking for negative values
            if self.__hour_representation__ == 'Normal':
                if self.__second__ < 0:
                    extra_min = self.__second__//-60
                    self.__minute__ -= 1 * extra_min+1
                    self.__second__ += 60 * (extra_min+1)
                if self.__minute__ < 0:
                    self.__hour__ -= 1
                    self.__minute__ += 60
            if self.__hour_representation__ != 'FOD':
                if self.__hour__ < 0:
                    self.__day__ -= 1
                    self.__hour__ += 24
            if self.__representation__ == 'DOY':
                if self.__day__ < 0:
                    self.__year__ -= 1
                    self.__day__ += 365 + leap_year(self.__year__)
            elif self.__representation__ == 'Gregorian':
                if self.__day__ < 1:
                    self.__month__ -= 1
                    if self.__month__ < 1:
                        self.__year__ -= 1
                        self.__month__ = 12
                    self.__day__ += days_in_month(self.__year__, self.__month__)

            # Start Checking for values that surpass the appropriate bounds
            if self.__hour_representation__ == 'Normal':
                if self.__second__ >= 60:
                    extra_min = self.__second__//60
                    self.__minute__ += extra_min
                    self.__second__ -= 60 * extra_min
                if self.__minute__ >= 60:
                    self.__hour__ += 1
                    self.__minute__ -= 60
            if self.__hour_representation__ != 'FOD':
                if self.__hour__ >= 24:
                    self.__day__ += 1
                    self.__hour__ -= 24
            if self.__representation__ == 'DOY':
                if self.__day__ > 365 + leap_year(self.__year__):
                    self.__year__ += 1
                    self.__day__ -= 365 + leap_year(self.__year__)
            elif self.__representation__ == 'Gregorian':
                if self.__day__ > days_in_month(self.__year__, self.__month__):
                    self.__day__ = self.__day__ - days_in_month(self.__year__, self.__month__)
                    self.__month__ += 1
                    if self.__month__ > 12:
                        self.__year__ += 1
                        self.__month__ = 1

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


# Updated Time class uses representation internally in modified Julian day, as to avoid boundary checking and then
# converts to the appropriate representation ONLY in the str method, to reduce the complexity of the class.
class Time:
    def __init__(self, time_sys, year=None, month=None, day=None, hour=None, minute=None, second=None):
        # Initialize all parameters of a time object
        self.__time_sys__ = time_sys
        self.__day__ = day

        # If the year is not inputted we assume the date is in Julian representation
        if year is None:
            self.__representation__ = 'Julian'
            self.__hour_representation__ = None

        # If the year is inputted and month is inputted we assume the date is in gregorian representation
        if year is not None and month is not None:
            self.__representation__ = 'Gregorian'
            self.__day__ = 367 * year - int((7 * (year + int((month + 9) / 12))) / 4) + int(
                275 * month / 9) + self.__day__ + 1721013.5

        # If the year is inputted and no month is inputted we assume the date is in date of year representation
        if year is not None and month is None:
            self.__representation__ = 'DOY'
            self.__day__ = 367 * year - int(7 * year / 4) + self.__day__ + 1721013.5

        # If the year is inputted and no hour is inputted we assume the date is in fraction of day representation
        if year is not None and hour is None:
            self.__hour_representation__ = 'FOD'

        # From this point on, the 'day' representation is assumed from the previous checks
        # If the hour is inputted and no minutes are inputted, we assume the time was inputted as a fraction of hour
        if year is not None and hour is not None and minute is None:
            self.__hour_representation__ = 'FOH'
            self.__day__ += hour/24
        # If a value is inputted for
        if year is not None and minute is not None:
            self.__hour_representation__ = 'Normal'
            self.__day__ += 1 / 24 * (hour + (1 / 60) * (minute + second / 60))

        if self.__day__ > 2400000.5:
            self.__day__ -= 2400000.5

    def __float__(self):
        return self.__day__

    # Overriding string method to allow proper formatting
    def __str__(self):
        if self.__representation__ == 'Julian':
            return str(self.__day__)

    def julian(self, modified=True):
        modified = not modified
        return self.__day__ + (2400000.5 * modified)

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

        # Add the time difference to the current time
        self.__day__ += time_diff / 86400
        self.__time_sys__ = 'TT'

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

        # Add the time difference to the current time
        self.__day__ += time_diff / 86400
        self.__time_sys__ = 'TAI'

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
        self.__day__ += time_diff / 86400
        self.__time_sys__ = 'UTC'

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

        # Add the time difference to the current time

        self.__day__ += time_diff / 86400
        self.__time_sys__ = 'UT1'

    def julian_cent(self):
        return (self.__day__ - 51544.5) / 36525

    def add_days(self, days):
        self.__day__ += days

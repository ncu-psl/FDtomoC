#include "common/time_process.h"

int isleap(int year) {
	int leap = 0;
	if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
		leap = 1;
	}
	return leap;
}

void etoh(double epoch, int *iyear, int *iday, int *ihour, int *imin,
	double *sec) {
	int dayinmon[13] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31 };
	/*
	 char monname[12][4] = { "Jan\0", "Feb\0", "Mar\0", "Apr\0", "May\0", "Jun\0",
	 "Jul\0", "Aug\0", "Sep\0", "Oct\0", "Nov\0", "Dec\0"
	 };
	 */
	int doy = (int)(epoch / 86400.0);
	double secleft = fmod(epoch, 86400.0);
	int hour = 0;
	int minute = 0;
	double second = 0.0;

	//---compute hours minutes seconds
	if (secleft != 0.0) {
		//---before 1970
		if (secleft < 0) {
			//---subtract a day
			doy = doy - 1;
			//----add a day
			secleft = secleft + 86400;
		}
		hour = (int)(secleft / 3600);
		secleft = fmod(secleft, 3600.0);
		minute = (int)(secleft / 60);
		second = fmod(secleft, 60.0);
	}

	int year, diy;
	if (doy >= 0) {
		year = 1970;
	a5: diy = 365 + isleap(year);
		if (doy < diy) {
			goto a10;
		}
		doy -= diy;
		year++;
		goto a5;
	}
	else {
		year = 1969;
	a7: diy = 365 + isleap(year);
		if (doy >= 0)
			goto a10;
		doy += diy;
		year--;
		goto a7;
	}
a10: doy++;
	//int date = year * 1000 + doy;
	int leap = isleap(year);
	int day = doy;
	int i = 0;
	for (i = 1; i < 12; i++) {
		int dim = dayinmon[i];
		if (leap && i == 1) {
			dim++;
		}
		if (day <= dim) {
			break;
		}
		day -= dim;
	}
	//int mon = i + 1;
	//har *mname = monname[i];

	*iyear = year;
	*iday = doy;
	*ihour = hour;
	*imin = minute;
	*sec = second;
}

void htoe(int iyear, int iday, int ihour, int imin, double sec, double *epoch) {
	int jdate = 1000 * iyear + iday;
	dtoepoch(jdate, epoch);
	*epoch = *epoch + 3600.0 * ihour + 60.0 * imin + sec;
}

void dtoepoch(int date, double *time) {
	int year = date / 1000;
	int day = date % 1000;
	int days = 0;
	if (year > 1970) {
		for (int i = 1970; i < year; i++) {
			days += 365;
			days += isleap(i);
		}
	}
	else if (year < 1970) {
		for (int i = year; i <= 1969; i++) {
			days -= 365;
			days -= isleap(i);
		}
	}
	days += (day - 1);
	*time = (double)(days * 24 * 60 * 60);
}

double htoe2(Time time){
	double epoch;
	int jdate = 1000 * time.iyr + time.jday;
	dtoepoch(jdate, &epoch);
	epoch = epoch + 3600.0 * time.ihr + 60.0 * time.imn + time.sec;
	return epoch;
}

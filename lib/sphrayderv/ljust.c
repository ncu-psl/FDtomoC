#include "sphrayderv/ljust.h"

//c-- - routine to left justify station names(no leading blanks)
//c-- - first make sure that string is not entirely blank
void ljust(char str[6]) {
	int i;
	for (i = 0; i < 6; i++)
		if (str[i] != ' ')
			break;
	while (str[0] == ' ') {
		for (i = 0; i < 5; i++) {
			str[i] = str[i + 1];
		}
		str[5] = ' ';
	}
}
/*
 *  utils.c — Clann v5.0.0
 *  Utility functions: numeric conversion, I/O helpers, printf2
 *
 *  Copyright (C) 2003-2026 Chris Creevey <chris.creevey@gmail.com>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "utils.h"

/* -----------------------------------------------------------------------
 * texttoint / inttotext
 * Simple single-character ↔ digit conversions used by toint / tofloat.
 * ----------------------------------------------------------------------- */

int texttoint(char c)
    {
    switch(c)
        {
        case '1' : return(1); break;
        case '2' : return(2); break;
        case '3':  return(3); break;
        case '4':  return(4); break;
        case '5':  return(5); break;
        case '6':  return(6); break;
        case '7':  return(7); break;
        case '8':  return(8); break;
        case '9':  return(9); break;
        default:   return(0); break;
        }
    }

char inttotext(int c)
    {
    switch(c)
        {
        case 1 : return('1'); break;
        case 2 : return('2'); break;
        case 3:  return('3'); break;
        case 4:  return('4'); break;
        case 5:  return('5'); break;
        case 6:  return('6'); break;
        case 7:  return('7'); break;
        case 8:  return('8'); break;
        case 9:  return('9'); break;
        default: return('0'); break;
        }
    }

/* -----------------------------------------------------------------------
 * toint / tofloat
 * Convert a decimal-digit string to int or float respectively.
 * Uses texttoint() for digit decoding and pow() from <math.h>.
 * ----------------------------------------------------------------------- */

int toint(char *number)
    {
    int charactercount = 0, j=0, result = 0;

    while(number[charactercount] != '\0') charactercount++;
    for(j=charactercount; j>0; j--)
        {
        result += (int)( pow(10, (charactercount-j)) * (texttoint(number[j-1])));
        }
    return(result);
    }

float tofloat(char *number)
    {
    int charactercount = 0, charactercount1 = 0, j=0;
    float real = 0, fraction = 0, total = 0;

    while(number[charactercount] != '.' && number[charactercount] != '\0') charactercount++;
    for(j=charactercount; j>0; j--)
        {
        real += (float)( pow(10, (charactercount-j)) * (texttoint(number[j-1])));
        }

    if(number[charactercount] == '.')
        {
        charactercount1 = charactercount;
        while(number[charactercount1] != '\0') charactercount1++;
        for(j=charactercount+1; j<charactercount1; j++)
            {
            fraction += (float)( pow(10,(charactercount-j))*(texttoint(number[j])));
            }
        }
    if(real == 0 && fraction == 0) total = 0;
    else total = real+fraction;
    return(total);
    }

/* -----------------------------------------------------------------------
 * xgets
 * A safe replacement for the deprecated gets() function.
 * Reads up to 79 characters from stdin, stopping at newline or CR.
 * Code pattern taken from "C: The Complete Reference" by Herbert Schildt.
 * ----------------------------------------------------------------------- */

char *xgets(char *s)
    {
    char ch, *p = NULL;
    int t = 0;

    for(t=0; t<80; t++)
        s[t] = '\0';
    p = s;

    for(t=0; t<80; ++t)
        {
        ch = getchar();
        switch(ch)
            {
            case '\n':
                s[t] = '\0';
                return(p);
                break;
            case '\r':
                s[t] = '\0';
                return(p);
                break;
            case '\b':
                if(t>0) t--;
                break;
            default:
                s[t] = ch;
                break;
            }
        }
    s[79] = '\0';
    return p;
    }

/* -----------------------------------------------------------------------
 * printf2
 * A drop-in replacement for printf that also writes to the log file when
 * logging is active (print_log == TRUE).  logfile and print_log are
 * defined in treecompare2.c and declared extern in utils.h.
 * ----------------------------------------------------------------------- */

void printf2(char *format, ...)
    {
    va_list ap;
    va_list ap2;

    if(print_log == TRUE)
        {
        va_start(ap, format);
        va_copy(ap2, ap);

        vfprintf(logfile, format, ap);
        va_end(ap);

        vprintf(format, ap2);
        va_end(ap2);
        }
    else
        {
        va_start(ap, format);
        vprintf(format, ap);
        va_end(ap);
        }
    }

/* A small C program that finds R from the registry, then runs a script */
/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena.c
 *
 * Description: This file provides a user interface to RSiena on Windows.
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
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses///
 * ****************************************************************************/

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>

static char *input[] ={
    "suppressPackageStartupMessages(library(RSiena))",
    "suppressPackageStartupMessages(library(tcltk))",
    "RSiena:::DONE(FALSE)",
    "siena01Gui()",
    "while(!RSiena:::DONE()) {Sys.sleep(0.01);}"};

#define REG_KEY_NAME "Software\\R-core\\R"

int main(int args, char **argv)
//int WINAPI WinMain(HINSTANCE hInst, HINSTANCE nn, LPSTR ll, int nCmdShow)
{
    char *p, tmp[1000], cmd[1000], errstr[1000];
    FILE *fp;
    int i, rc;
    DWORD type, size = 1000;
    BYTE data[1000], *d = data;
    HKEY hkey;
    PROCESS_INFORMATION pi;
    STARTUPINFO si;
    int ret;

    p = getenv("TEMP");
    if(!p) p = getenv("TMP");
    if(!p) {
		sprintf(errstr, "neither TEMP nor TMP is set\n");
		MessageBox(NULL, errstr, NULL, 0);
		exit(2);
    }
	sprintf(tmp, "%s/RSiena.tmp", p);

	  fp = fopen(tmp, "wt");
    if (!fp) {
		sprintf(errstr, "cannot open temporary file %s\n", tmp);
		MessageBox(NULL, errstr, NULL, 0);
		exit(3);
    }
    for (i = 0; i < sizeof(input)/sizeof(char *); i++)
		fprintf(fp,"%s\n", input[i]);
    fclose(fp);

    if ((rc = RegOpenKeyEx(HKEY_LOCAL_MACHINE, REG_KEY_NAME, 0,
				KEY_READ, &hkey)) != ERROR_SUCCESS) {
        if ((rc = RegOpenKeyEx(HKEY_CURRENT_USER, REG_KEY_NAME, 0,
					KEY_READ, &hkey)) != ERROR_SUCCESS) {

            sprintf(errstr,"R is not registered\n");
			MessageBox(NULL, errstr, NULL, 0);
				exit(1);
		}}
    rc = RegQueryValueEx(hkey, "installPath", NULL, &type, d, &size);
    if (rc != ERROR_SUCCESS) {
		sprintf(errstr,"No current version of R is registered\n");
		MessageBox(NULL, errstr, NULL, 0);
		exit(1);
    }

    sprintf(cmd, "%s\\bin\\Rterm.exe --slave --no-init-file -f %s",
(char *)d, tmp);
    printf("cmd is %s\n", cmd);
    ZeroMemory( &si, sizeof(si) );
    ZeroMemory( &pi, sizeof(pi) );

    si.cb = sizeof(si);
    si.lpReserved = NULL;
    si.lpReserved2 = NULL;
    si.cbReserved2 = 0;
    si.lpDesktop = NULL;
    si.lpTitle = NULL;
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_HIDE;
    ret = CreateProcess(NULL, cmd, NULL, NULL, FALSE,
                        0, NULL, NULL, &si, &pi);
    exit(0);
}

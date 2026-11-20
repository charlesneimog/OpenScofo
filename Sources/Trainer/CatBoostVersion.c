/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.

    CatBoost's ONNX exporter records these build-information symbols.
*/

#include <stdio.h>

#ifndef OPENSCOFO_CATBOOST_VERSION
#define OPENSCOFO_CATBOOST_VERSION "unknown"
#endif

const char *GetProgramSvnVersion(void) {
    return "CatBoost " OPENSCOFO_CATBOOST_VERSION;
}

void PrintProgramSvnVersion(void) {
    puts(GetProgramSvnVersion());
}

const char *GetArcadiaSourcePath(void) {
    return "";
}

const char *GetArcadiaSourceUrl(void) {
    return "";
}

const char *GetArcadiaLastChange(void) {
    return "";
}

int GetArcadiaLastChangeNum(void) {
    return 0;
}

const char *GetArcadiaLastAuthor(void) {
    return "";
}

int GetProgramSvnRevision(void) {
    return 0;
}

const char *GetVCSDirty(void) {
    return "";
}

const char *GetProgramHash(void) {
    return "";
}

const char *GetProgramCommitId(void) {
    return "";
}

const char *GetCustomVersion(void) {
    return "";
}

const char *GetReleaseVersion(void) {
    return OPENSCOFO_CATBOOST_VERSION;
}

const char *GetProgramScmData(void) {
    return "";
}

const char *GetProgramShortVersionData(void) {
    return OPENSCOFO_CATBOOST_VERSION;
}

const char *GetProgramBuildUser(void) {
    return "";
}

const char *GetProgramBuildHost(void) {
    return "";
}

const char *GetProgramBuildDate(void) {
    return "";
}

int GetProgramBuildTimestamp(void) {
    return 0;
}

const char *GetVCS(void) {
    return "git";
}

const char *GetBranch(void) {
    return "";
}

const char *GetTag(void) {
    return OPENSCOFO_CATBOOST_VERSION;
}

int GetArcadiaPatchNumber(void) {
    return 0;
}

void PrintSvnVersionAndExit(int argc, char *argv[]) {
    (void)argc;
    (void)argv;
}

void PrintSvnVersionAndExitEx(int argc, char *argv[], const char *opts) {
    (void)argc;
    (void)argv;
    (void)opts;
}

void PrintSvnVersionAndExit0(void) {
}

/*

    MVPTree c library
    Copyright (C) 2008-2009 Aetilius, Inc.
    All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    D Grant Starkweather - dstarkweather@phash.org

*/

// imget.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <windows.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <pHash.h>

#define strncasecmp _strnicmp

extern "C" {
#include "mvptree.h"
}

#define MVP_BRANCHFACTOR 2
#define MVP_PATHLENGTH   5
#define MVP_LEAFCAP     25

typedef int(*HashFunc)(const char* file, ulong64 &hash);

static uint64_t nbcalcs = 0;

float hamming_distance(MVPDP* pointA, MVPDP* pointB)
{
    if(!pointA || !pointB || pointA->datalen != pointB->datalen) return -1.0f;

    uint64_t a = *((uint64_t*)pointA->data);
    uint64_t b = *((uint64_t*)pointB->data);
    uint64_t x = a ^ b;
    
    const uint64_t m1  = 0x5555555555555555ULL;
    const uint64_t m2  = 0x3333333333333333ULL;
    const uint64_t h01 = 0x0101010101010101ULL;
    const uint64_t m4  = 0x0f0f0f0f0f0f0f0fULL;

    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;

    float result = (float)((x * h01) >> 56);
    result = exp(result - 1);
    nbcalcs++;
    return result;
}


int main(int argc, char* argv[])
{
    if(argc < 3) {

        printf("not enough input args\n");
        printf("%Fs command filename [directory] [radius] [k-nearest]\n\n", argv[0]);
        printf("  command    - command - e.g. 'add', 'query' or 'print'\n");
        printf("  filename   - file from which to read the tree\n");
        printf("  directory  - directory (for add and query)\n");
        printf("  radius     - radius for query operation (default = 21.0)\n");
        printf("  k-nearest  - max number of results (default = 5)\n");
        return 0;
    }

    const char* command  = argv[1];
    const char* filename = argv[2];
    const char* dirname  = (argc > 3) ? argv[3] : "None";
    const float radius   = (argc > 4) ? atof(argv[4]) : 21.0;
    const int knearest = (argc > 5) ? atoi(argv[5]) : 5;
    const int blur = (argc > 6) ? atoi(argv[6]) : NULL;
    
    printf("command   - %s\
            \nfilename  - %s\
            \ndir       - %s\
            \nradius    - %g\
            \nk-nearest - %d\n",
            command, filename, dirname, radius, knearest);

    CmpFunc distance_func = hamming_distance;
    int nbfiles = 0;
    char** files = NULL;
    MVPDP** points;

    MVPError err;
    MVPTree* tree = mvptree_read(filename, distance_func, MVP_BRANCHFACTOR, MVP_PATHLENGTH, MVP_LEAFCAP, &err);
    assert(tree);

    if(!strncasecmp(command, "add", 3) || !strncasecmp(command, "query", 3)) {

        if((_access( dirname, 0 )) == -1)
        {
            perror("\nError" );
            return 0;
        }

        files = ph_readfilenames(dirname, nbfiles);
        if(nbfiles == 0){
            printf("\nError: No valid files in %s\n", dirname);
            return 0;
        }
        printf("\n%d files in %s\n\n", nbfiles, dirname);

        points = (MVPDP**)malloc(nbfiles * sizeof(MVPDP*));

        int count = 0;
        uint64_t hashvalue;

        clock_t begin, end;
        double time_spent;

        HashFunc hashfunction = (blur == NULL) ? ph_dct_imagehash_no_blur : ph_dct_imagehash;

        begin = clock();
        for(int i = 0; i < nbfiles; i++) {
            if(hashfunction(files[i], hashvalue) < 0) {
                printf("Unable to get hash value for %s.\n", files[i]);
                continue;
            }

            printf("%016llx %s\n", hashvalue, files[i]);

            points[count] = dp_alloc(MVP_UINT64ARRAY);
            points[count]->id = strdup(files[i]);
            points[count]->data = malloc(1 * MVP_UINT64ARRAY);
            // our data (hash) will be 1 MVP_UINT64ARRAY (8 bytes) in size
            points[count]->datalen = 1;
            memcpy(points[count]->data, &hashvalue, MVP_UINT64ARRAY);
            count++;
        }

        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

        printf("\n");
        printf("Time: %.2f\n", time_spent);

        if(!strncasecmp(command, "add", 3)) {
            //printf("\nAdding %d hashes to tree.\n", count);

            MVPError error = mvptree_add(tree, points, count);

            if(error != MVP_SUCCESS) {
                printf("Unable to add hash values to tree: %s\n", mvp_errstr(error));
                goto cleanup;
            }

            error = mvptree_write(tree, filename, S_IWRITE);

            if(error != MVP_SUCCESS) {
                perror("Unable to save file");
                goto cleanup;
            }
        }

        else if(!strncasecmp(command, "query", 3)) {
            unsigned int nbresults = 0, j = 0;
            FILE *stream;

            fopen_s( &stream, "query.html", "w" );

            printf("Writing query.html\n");

            fprintf(stream, "<!doctype html> \
                    <html lang=\"en\"> \
                     <head> \
                      <meta charset=\"UTF-8\"> \
                      <title> Image Query </title> \
                      <style> \
                      .q {padding: auto 15px} \
                      img {max-width: 17%; max-height: 250px;} \
                      div {padding: 15px; border-bottom: 1px dashed black;} \
                      </style> \
                     </head> \
                     <body>\n");

            for(int i = 0; i < count; i++) {
                nbcalcs = 0;
                MVPDP** results = mvptree_retrieve(tree, points[i], knearest, radius, &nbresults, &err);

                fprintf(stream, "<div><img src='%s' title='%s' class='q'>", files[i], files[i]);

                if(nbresults > 0) {
                    for(j = 0; j < nbresults; j++) {
                        fprintf(stream, "\n<img src='%s' title='%s'>", results[j]->id, results[j]->id);
                    }
                } else {
                        fprintf(stream, "\nNo matching images found.");
                    }
                

                fprintf(stream, "</div>\n");

                free(results);
            }

            fprintf(stream, " </body></html>");

        }
    }

    else if(!strncasecmp(command, "print", 3)) {
        printf("-----------------------print-------------------------\n");
        mvptree_print(stdout, tree);
        printf("-----------------------------------------------------\n\n");
    }

    mvptree_clear(tree, free);

cleanup:

    if(nbfiles) {
        for(int i = 0; i < nbfiles; i++) {
            free(files[i]);
        }
    }

    if(files) free(files);

    return 0;
}

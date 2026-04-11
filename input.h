#ifndef INPUT_H
#define INPUT_H

typedef struct INPUT
{
    int N;
    int Nmax;
    char** name;
    char** value;

} INPUT;

INPUT* inputInit(char* fileName, int N);

void inputPrint(INPUT* input);

int inputNameIsInput(INPUT* input, char* name);

char* inputGetValue(INPUT* input, char* name);

void inputFree(INPUT* input);

#endif

#ifndef INPUT_H
#define INPUT_H

typedef struct COMMAND
{

    char* variable;
    char* value;
    struct COMMAND* next;
    
} COMMAND;


typedef struct INPUT
{

    COMMAND* firstCommand;
    int Ncommand;

} INPUT;


INPUT* inputInit(char* fileName);

COMMAND* inputCommandInit(char* variable, char* value);

void inputCommandFree(COMMAND* command);

void inputCommandPrint(COMMAND* command);

void inputPrint(INPUT* input);

int inputNameIsInput(INPUT* input, char* name);

char* inputGetValue(INPUT* input, char* name);

void inputFree(INPUT* input);

#endif

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include"input.h"

INPUT* inputInit(char* fileName)
{

    INPUT* input = malloc(sizeof(INPUT));
    FILE* ff = fopen(fileName, "r");
    char c;
    
    char variable[50];
    char value[50];
    
    variable[0] = '\0';
    value[0] = '\0';
    
    bool getVar = true;
    bool getValue = false;
    bool save = true;
    bool first = true;
    
    COMMAND* command = NULL;

    int nVar = 0;
    int nValue = 0;

    input->Ncommand = 0;
    
    c = fgetc(ff);    
    
    while(c != EOF)
    {   
        if(c == ' ' || c == ',' || c == '\n' || c == '\r' || c == '#')
        {
            save = false;
        }
        else
        {
            save = true;
        }
        
        if(c == '\n' || c == '\r')
        {
            if(nVar > 0 && nValue > 0)
            {
                if(first)
                {
                    command = inputCommandInit(variable, value);
                    input->firstCommand = command;
                    first = false;
                }
                else
                {
                    command->next = inputCommandInit(variable, value);
                    command = command->next;
                }
                
                input->Ncommand += 1;
            }
                        
            getVar = true;
            getValue = false;
            
            variable[0] = '\0';
            value[0] = '\0';   
            
            nVar = 0;
            nValue = 0;
        }

        if(c == ',')
        {
            getVar = false;
            getValue = true;
        }

        if(c == '#')
        {
            getVar = false;
            getValue = false;
        }
     
        if(getVar)
        {
            if(save)
            {                  
                variable[nVar] = c;
                variable[nVar + 1] = '\0';
                nVar += 1;
            }
        }
         
        if(getValue)
        {
            if(save)
            {            
                value[nValue] = c;
                value[nValue + 1] = '\0';
                nValue += 1;
            }
        }

        c = fgetc(ff);
    }
    
    fclose(ff);
    
    return input;
}


COMMAND* inputCommandInit(char* variable, char* value)
{
    COMMAND* command = malloc(sizeof(COMMAND));
    command->variable = malloc(50*sizeof(char));
    command->value = malloc(50*sizeof(char));

    strcpy(command->variable, variable);
    strcpy(command->value, value);    
    
    command->next = NULL;
        
    return command;
}


void inputCommandFree(COMMAND* command)
{
    free(command->variable);
    free(command->value);
    free(command);
}

void inputCommandPrint(COMMAND* command)
{
    printf(" %s, %s\n", command->variable, command->value);
}


void inputPrint(INPUT* input)
{
    COMMAND* command = input->firstCommand;
    while(command)
    {
        inputCommandPrint(command);
        command = command->next;
    }
}


int inputNameIsInput(INPUT* input, char* name)
{
    
    int found = 0;
    
    COMMAND* command = input->firstCommand;
    while(command)
    {               
        if(strcmp(command->variable, name) == 0)
        {
            found = 1;
            break;
        }
        command = command->next;
    }
    
    return found;
    
}


char* inputGetValue(INPUT* input, char* name)
{
    char* s;
    int found = 0;
    
    COMMAND* command = input->firstCommand;
    while(command)
    {               
        if(strcmp(command->variable, name) == 0)
        {
            s = command->value;
            found = 1;
            break;
        }
        command = command->next;
    }
    
    if(found == 0)
    {
        printf("Error: Value not found in the input file: %s.\n", name);
        exit(0);
    }

    return s;
    
}


void inputFree(INPUT* input)
{
    COMMAND* command = input->firstCommand;
    COMMAND* next = command->next;
    
    while(next)
    {
        inputCommandFree(command);
        command = next;
        next = command->next;
    }
    
    inputCommandFree(command);
    
    free(input);    
} 


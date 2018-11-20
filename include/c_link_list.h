#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LINKED_LIST_NULL -1
#define NOT_FOUND -2

typedef struct NODES {
	char *key;
	struct NODES *next;
} Node;

Node * create_new_node(char *key);
void link_list_add(Node *, Node *);
int link_list_get_index(Node *, char *);
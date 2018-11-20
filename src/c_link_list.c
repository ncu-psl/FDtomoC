#include "c_link_list.h"

Node * create_new_node(char *key) {
	Node *new_node = (Node *)malloc(sizeof(Node));
	new_node->key=key;
	new_node->next=NULL;
	return new_node;
}

void link_list_add(Node *list_head, Node *node) {
	Node *ptr=list_head;
	while(ptr->next!=NULL) {
		ptr=ptr->next;
	}
	ptr->next = node;
}

int link_list_get_index(Node *list_head, char *str) {
	if(list_head==NULL) {
		return LINKED_LIST_NULL;
	}

	int index = 0;
	Node *ptr=list_head;
	while(ptr) {
		if(strcmp(ptr->key,str)==0) {
			return index;
		}
		ptr=ptr->next;
		index++;
	}
	return NOT_FOUND;
}

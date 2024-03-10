#include <stdio.h>
// ---------------------------------------------------------------------------
#include "list_void_ptr.h"
// ---------------------------------------------------------------------------
TLista::TLista()
{
	head = NULL;
	tail = NULL;
	count = 0;
}
// ---------------------------------------------------------------------------
TLista::TLista(TLista& c)
{
	head = NULL;
	tail = NULL;

	node_double_linked *tmp = c.head;
	while (tmp)
	{
		Add(tmp->inf);
		tmp = tmp->next;
	}
	count = c.count;
}
// ---------------------------------------------------------------------------
TLista::~TLista()
{
	while (head)
	{
		node_double_linked * tmp = head;
		head = head->next;
		delete tmp;
	}
	tail = NULL;
	count = 0;
}

// ---------------------------------------------------------------------------
void TLista::Add(void* data)
{
	if (!head){
		head = new node_double_linked;
		head->inf = data;
		head->next = NULL;
		head->prev = NULL;
		tail = head;
	}
	else
	{
		node_double_linked* tmp = new node_double_linked;
		tmp->inf = data;
		tmp->next = NULL;
		tmp->prev = tail;
		tail->next = tmp;
		tail = tmp;
	}
	count++;
}
// ---------------------------------------------------------------------------
void TLista::Delete(int Index)
{
	if (Index < count)
		if (head)
		{
			if (!Index)
			{ // delete the head
				node_double_linked* tmp = head;
				head = head->next;
				if (head)
					head->prev = NULL;
				delete tmp;
				if (count == 1)
					tail = NULL;
			}
			else
			{ // delete the other
				int i = 0;
				node_double_linked* tmp = head;
				while (tmp && (i < Index - 1))
				{
					i++;
					tmp = tmp->next;
				}
				// sunt pozitionat pe anterioru
				if (tmp)
				{
					node_double_linked* urm = tmp->next;
					tmp->next = urm->next;
					if (urm->next)
						urm->next->prev = tmp;
					delete urm;
					if (Index == count - 1)
						tail = tmp;
				}
			}
			count--;
		}
}

// ---------------------------------------------------------------------------
void* TLista::GetInfo(int Index)
{
	if (!head)
		return NULL;

	int i = 0;
	node_double_linked* tmp = head;
	while (tmp && (i < Index))
	{
		i++;
		tmp = tmp->next;
	}
	if (tmp)
		return tmp->inf;
	else
		return NULL;
}

// ---------------------------------------------------------------------------
node_double_linked* TLista::GetNode(int Index)
{
	if (!head)
		return NULL;

	int i = 0;
	node_double_linked* tmp = head;
	while (tmp && (i < Index))
	{
		i++;
		tmp = tmp->next;
	}
	if (tmp)
		return tmp;
	else
		return NULL;
}

// ---------------------------------------------------------------------------
void* TLista::GetNextInfo(node_double_linked* p)
{
	if (p){
		if (p->next)
			return p->next->inf;
		else
			return NULL;
    }
	return NULL;
}

// ---------------------------------------------------------------------------
void* TLista::GetPrevInfo(node_double_linked* p)
{
	if (p){
		if (p->prev)
			return p->prev->inf;
		else
			return NULL;
    }
	return NULL;
}
// ---------------------------------------------------------------------------
void* TLista::GetNextCircularInfo(node_double_linked* p)
{
	if (p){
		if (!p->next)
			return head->inf;
		else
			return p->next->inf;
    }
	return NULL;
}

// ---------------------------------------------------------------------------
void* TLista::GetPrevCircularInfo(node_double_linked* p)
{
	if (p){
		if (!p->prev)
			return tail->inf;
		else
			return p->prev->inf;
    }
	return NULL;
}

// ---------------------------------------------------------------------------
node_double_linked* TLista::delete_current_circular(node_double_linked* p)
{
	if (count > 1) {
		// more elements
		if (p == head)
			head = head->next;
		if (p == tail)
			tail = tail->prev;

		p->prev->next = p->next;
		p->next->prev = p->prev;
		delete p;
		count--;
		return head;
	}
	else {
		// one element
		delete p;
		tail = NULL;
		head = NULL;
		count = 0;
		return NULL;
	}
}
// ---------------------------------------------------------------------------
node_double_linked* TLista::DeleteCurrent(node_double_linked* p)
{
	if (count > 0)
		if (p == head)
		{
			head = head->next;
			if (head)
				head->prev = NULL;
			delete p;
			// p = NULL;
			if (count == 1)
				tail = NULL;
			count--;
			return head;
		}
		else
			if (p == tail)
			{
				tail = tail->prev;
				if (tail)
					tail->next = NULL;
				delete p;
				// p = NULL;
				if (count == 1)
					head = NULL;
				count--;
				return NULL;
			}
			else
			{
				node_double_linked *tmp = p->next;
				p->prev->next = p->next;
				p->next->prev = p->prev;
				delete p;
				count--;
				return tmp;
			}
	else
		return NULL;
}

// ---------------------------------------------------------------------------
void TLista::Append(TLista & )
{

}

// ---------------------------------------------------------------------------
void TLista::Insert(long Index, void* data)
{

	if (Index >= count)
		Add(data);
	else
		if (Index == 0)
		{ // insert in the head
			node_double_linked* k = new node_double_linked;
			k->inf = data;
			k->prev = NULL;
			k->next = head;
			head->prev = k;
			head = k;
			count++;
		}
		else
		{
			long i = 0;
			node_double_linked* tmp = head;
			while (tmp && (i < Index - 1))
			{
				i++;
				tmp = tmp->next;
			}
			// insert AFTER tmp
			if (tmp)
			{
				node_double_linked* k = new node_double_linked;
				k->inf = data;
				k->prev = tmp;
				k->next = tmp->next;
				tmp->next->prev = k;
				tmp->next = k;
				count++;
			}
		}
}

// ---------------------------------------------------------------------------
void TLista::Clear(void)
{
	while (head)
	{
		node_double_linked * tmp = head;
		head = head->next;
		delete tmp;
	}
	tail = NULL;
	head = NULL;
	count = 0;
}

// ---------------------------------------------------------------------------
void TLista::DeleteHead(void)
{
	if (head)
	{
		node_double_linked* tmp = head;
		head = head->next;
		if (head)
			head->prev = NULL;
		delete tmp;
		if (count == 1)
			tail = NULL;
		count--;
	}
}

// ---------------------------------------------------------------------------
	void* TLista::GetCurrentInfo(node_double_linked* p){
		if (p)
			return p->inf;
		return 0;
	}
// ---------------------------------------------------------------------------
	void* TLista::GetHeadInfo(void){
		if (head)
			return head->inf;
		return 0;
	}
// ---------------------------------------------------------------------------
	void* TLista::GetTailInfo(void)
	{
		if (tail)
			return tail->inf;
		return 0;
	}
	// ---------------------------------------------------------------------------
	void AppendWithoutCopy(TLista &)
{

}
// ---------------------------------------------------------------------------
void TLista::make_circular(void)
{
	if (head) {
		head->prev = tail;
		tail->next = head;
	}
}
// ---------------------------------------------------------------------------

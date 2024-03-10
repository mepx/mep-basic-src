#ifndef lista_void_ptr_H
#define lista_void_ptr_H
//---------------------------------------------------------------------------
// simple linked dynamically allocated list
//---------------------------------------------------------------------------
struct node_double_linked{
	void* inf;
	node_double_linked* next, *prev;
};
//---------------------------------------------------------------------------
class TLista{
private:
public:
	node_double_linked* head, *tail;
	int count;

	TLista();
	TLista(TLista&);
	~TLista();
	void Add(void*);
	void Delete(int);
	node_double_linked* DeleteCurrent(node_double_linked*);
	void* GetInfo(int);
	node_double_linked* GetNode(int Index);
	void Append(TLista & source);
	void Insert(long, void*);
	void Clear(void);
	void DeleteHead(void);

	void* GetCurrentInfo(node_double_linked* p);
	void* GetHeadInfo(void);
	void* GetTailInfo(void);

	void* GetNextInfo(node_double_linked*);
	void* GetPrevInfo(node_double_linked*);
};
//---------------------------------------------------------------------------
#endif

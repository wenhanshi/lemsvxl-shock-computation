// This is brcv/shp/dbsk2d/dbsk2d_ishock_graph.cxx

//:
// \file

#include "dbsk2d_ishock_graph.h"

//: simpler constructor by Wenhan
dbsk2d_ishock_graph::dbsk2d_ishock_graph()
{
    _nextAvailableID = 0;
    _boundary = NULL;
    edges_.clear();
    vertices_.clear();
}

//: Constructor
dbsk2d_ishock_graph::dbsk2d_ishock_graph(dbsk2d_boundary_sptr boundary)
{
  _nextAvailableID = 0;
  _boundary = boundary;

  edges_.clear();
  vertices_.clear();

  if (_boundary){
    // also link the boundary object to this shock objact
    _boundary->set_ishock_graph(this);
  }
}

//: Destructor
dbsk2d_ishock_graph::~dbsk2d_ishock_graph ()
{
  if (_boundary)
    _boundary->set_ishock_graph(NULL);

  clear();
}

//-------------------------------------------------------------------
// Access operations
//-------------------------------------------------------------------

//: returns the number of nodes in shock graph.
int dbsk2d_ishock_graph::number_of_nodes() const
{
  return vertices_.size();
}

//: returns the number of edges in shock graph.
int dbsk2d_ishock_graph::number_of_edges() const
{ 
  return edges_.size();
}

//: returns the number of edges adjacent to node v (| adj_edges(v)|).
int dbsk2d_ishock_graph::outdeg(dbsk2d_ishock_node* v, bool exclude_hidden)
{
  return v->outdeg(exclude_hidden);
}

//: returns the number of edges ending at v (| in_edges(v)|)
int dbsk2d_ishock_graph::indeg(dbsk2d_ishock_node* v, bool exclude_hidden)
{
  return v->indeg(exclude_hidden);
}

//: returns outdeg(v) + indeg(v).
int dbsk2d_ishock_graph::degree(dbsk2d_ishock_node* v, bool exclude_hidden)
{
  return v->degree(exclude_hidden);
}

//: returns adj_edges(v). 
// \todo {finish this.}
ishock_edge_list dbsk2d_ishock_graph::adj_edges(dbsk2d_ishock_node* v)
{
  ishock_edge_list linksList;
  return linksList;
}

//: returns the list of all nodes adjacent to v.
// \todo {finish this.} 
ishock_node_list dbsk2d_ishock_graph::adj_nodes(dbsk2d_ishock_node* v)
{
  ishock_node_list nodesList;
  return nodesList;
}

//-------------------------------------------------------------------
// Graph building functions
//-------------------------------------------------------------------

//: Adds a new vertex to the graph
void dbsk2d_ishock_graph::add_vertex(dbsk2d_ishock_node* vertex)
{
  if (!vertex) return;

  vertices_.push_back(vertex);
}

//: Deletes a vertex in the graph
void dbsk2d_ishock_graph::remove_vertex(dbsk2d_ishock_node* vertex)
{
  if (!vertex) return;

  vertices_.remove(vertex);

  //since these are not smart pointers we have to actively delete them
  delete vertex;
}

//: Add an edge
void dbsk2d_ishock_graph::add_edge( dbsk2d_ishock_edge* edge)
{
  if (!edge) return;

  edges_.push_back(edge);

  //Add this shock to both the boundary elements
  //(i.e., make new wavefront branches on the wavefront trees)
  edge->lBElement()->add_shock (edge);
  edge->rBElement()->add_shock (edge);
}


//: Remove an edge
void dbsk2d_ishock_graph::remove_edge( dbsk2d_ishock_edge* edge)
{
  if (!edge) return;

  edges_.remove(edge);

  //delete the links from the boundary elements to this edge
  bool success = edge->lBElement()->delete_shock (edge) &&
                 edge->rBElement()->delete_shock (edge);
  dbsk2d_assert(success);

  if ( edge->rNeighbor())
  {
      if ( edge->rNeighbor()->lShock() == edge )
      {
          edge->rNeighbor()->set_lShock(NULL);
      }
      else if( edge->rNeighbor()->rShock() == edge )
      {
          edge->rNeighbor()->set_rShock(NULL);
      }
      
      edge->rNeighbor()->clear_lNeighbor();
  }

  if ( edge->lNeighbor())
  {
      if ( edge->lNeighbor()->rShock() == edge )
      {
          edge->lNeighbor()->set_rShock(NULL);
      }
      else if( edge->lNeighbor()->lShock() == edge )
      {
          edge->lNeighbor()->set_lShock(NULL);
      }

      edge->lNeighbor()->clear_rNeighbor();
  }

  edge->lBElement()->delete_lr_refs(edge);
  edge->rBElement()->delete_lr_refs(edge);
  
  //since these are not smart pointers we have to actively delete them
  delete edge;
  edge=0;
}

void dbsk2d_ishock_graph::clear()
{
  //delete all the nodes and edges
  while(edges_.size()>0)
    remove_edge(edges_.back());

  while(vertices_.size()>0)
    remove_vertex(vertices_.back());

  //reset the id counter
  _nextAvailableID = 0;
}

//-------------------------------------------------------------------
// Directed graph operations
//-------------------------------------------------------------------

//: returns the source node of edge e.
dbsk2d_ishock_node* dbsk2d_ishock_graph::source(dbsk2d_ishock_edge* e)
{
  return e->pSNode();
}

//: returns the dest node of edge e.
dbsk2d_ishock_node* dbsk2d_ishock_graph::dest(dbsk2d_ishock_edge* e)
{
  return e->cSNode();
}

//:  returns dest(e) if v = source(e) and source(e) if v = dest(e)
// if either of these are NULL it returns itself
// because the assumption is that we are turning around 
// a semi infinite edge of the graph
dbsk2d_ishock_node* dbsk2d_ishock_graph::opposite(dbsk2d_ishock_node* v, dbsk2d_ishock_edge* e)
{
  
  dbsk2d_ishock_node* oppositeNode;

  if (v==e->pSNode())
    oppositeNode = e->cSNode();
  else
    oppositeNode = e->pSNode();

  if (!oppositeNode)
    return v;
  else
    return oppositeNode;
}

//:  returns adj_edges(v) if dbsk2d_ishock_graph is directed and the empty list otherwise.
// \todo {finish this.}
ishock_edge_list dbsk2d_ishock_graph::out_edges(dbsk2d_ishock_node* v)
{
  ishock_edge_list linksList;
  return linksList;
}

//:  returns in_edges(v) if dbsk2d_ishock_graph is directed and the empty list otherwise.
// \todo {finish this.}
ishock_edge_list dbsk2d_ishock_graph::in_edges(dbsk2d_ishock_node* v)
{
  ishock_edge_list linksList;
  return linksList;
}

//:  returns the first edge in the adjacency list of v (nil if this list is empty). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::first_adj_edge(dbsk2d_ishock_node* v)
{
  return NULL;
}

//:  returns the last edge in the adjacency list of v (nil if this list is empty).
// \todo {finish this.} 
dbsk2d_ishock_edge* dbsk2d_ishock_graph::last_adj_edge(dbsk2d_ishock_node* v)
{
  return NULL;
}

//:  returns the successor of edge e in the adjacency list of node source(e) (nil if it does not exist). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::adj_succ(dbsk2d_ishock_edge* e)
{
  return NULL;
}

//:  returns the predecessor of edge e in the adjacency list of node source(e) (nil if it does not exist). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::adj_pred(dbsk2d_ishock_edge* e)
{
  return NULL;
}

//:  returns the cyclic successor of edge e in the adjacency list of node source(e). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::cyclic_adj_succ(dbsk2d_ishock_edge* e)
{
  return NULL;
}

//:  returns the cyclic predecessor of edge e in the adjacency list of node source(e). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::cyclic_adj_pred(dbsk2d_ishock_edge* e)
{
  return NULL;
}

//v returns the first edge of in_edges(v) (nil if this list is empty). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::first_in_edge(dbsk2d_ishock_node* v)
{
  return NULL;
}

//:  returns the last edge of in_edges(v) (nil if this list is empty). 
// \todo {finish this.} 
dbsk2d_ishock_edge* dbsk2d_ishock_graph::last_in_edge(dbsk2d_ishock_node* v)
{
  return NULL;
}

//:  returns the successor of edge e in in_edges(target(e)) (nil if it does not exist). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::in_succ(dbsk2d_ishock_edge* e)
{
  return NULL;
}

//:  returns the predecessor of edge e in in_edges(target(e)) (nil if it does not exist).
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::in_pred(dbsk2d_ishock_edge* e)
{
  return NULL;
}

//:  returns the cyclic successor of edge e in in_edges(target(e)) (nil if it does not exist). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::cyclic_in_succ(dbsk2d_ishock_edge* e)
{
  return NULL;
}

//:  returns the cyclic predecessor of edge e in in_edges(target(e)) (nil if it does not exist). 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::cyclic_in_pred(dbsk2d_ishock_edge* e)
{
  return NULL;
}

//-------------------------------------------------------------------
// Undirected graph operations   
//-------------------------------------------------------------------

//:  returns the successor of edge e in the adjacency list of v.
//Precondition: e is incident to v. 
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::adj_succ(dbsk2d_ishock_edge* e, dbsk2d_ishock_node* v)
{
  return NULL;
}

//:  returns the predecessor of edge e in the adjacency list of v.
//Precondition: e is incident to v.
// \todo {finish this.}
dbsk2d_ishock_edge* dbsk2d_ishock_graph::adj_pred(dbsk2d_ishock_edge* e, dbsk2d_ishock_node* v)
{  
  return NULL;
}

//:  returns the cyclic successor(CW) of edge e in the adjacency list of v.
// Precondition: e is incident to v. 
// if this is a leaf edge, it should return itself
dbsk2d_ishock_edge* dbsk2d_ishock_graph::cyclic_adj_succ(dbsk2d_ishock_edge* e, dbsk2d_ishock_node* v, bool exclude_hidden)
{
  if (v->is_a_source()){
    //assume that a source node cannot exist with just one child
    if (e == ((dbsk2d_ishock_node*)v)->cShock())
      return ((dbsk2d_ishock_node*)v)->cShock2();
    else
      return ((dbsk2d_ishock_node*)v)->cShock();
  }
  else if (v->is_an_A3source()){
    //an A3 is always a leaf edge
    return e;
  }
  else {
    //junctions and sinks can have multiple parents 
    //start looking for the current edge starting from child
    //any one of the edges incident on these nodes can be hidden

    bool bEdgeFound = false;

    dbsk2d_ishock_edge* cSLink = NULL;
    if (v->is_a_junct())
      if ( ((dbsk2d_ishock_node*)v)->cShock())
        if ( ((dbsk2d_ishock_node*)v)->cShock()->isNotHidden() || !exclude_hidden )
          cSLink = ((dbsk2d_ishock_node*)v)->cShock();

    if (e == cSLink)
      bEdgeFound = true;

    //now go through the list of parents
    ishock_edge_list::iterator curS = v->pShocks().begin();
    for(; curS!=v->pShocks().end(); curS++){
      dbsk2d_ishock_edge* curSLink = (*curS);
      if (curSLink->isHidden() && exclude_hidden) curSLink = NULL;

      if (bEdgeFound && curSLink)
        return curSLink;

      if (e == curSLink)
        bEdgeFound = true;
    }

    if (bEdgeFound){ //this should always be true by now
      if (cSLink)
        return cSLink;
      else {
        //need to go around the parents one more time
        ishock_edge_list::iterator curS = v->pShocks().begin();
        for(; curS!=v->pShocks().end(); ++curS){
          dbsk2d_ishock_edge* curSLink = (*curS);
          if (curSLink->isNotHidden() || !exclude_hidden) 
            return curSLink;
        }
      }
    }
  }
  
  //it should never get to here
  return NULL;
} 

//:  returns the cyclic predecessor(CW) of edge e in the adjacency list of v.
// Precondition: e is incident to v. 
// if this is a leaf edge, it should return itself
dbsk2d_ishock_edge* dbsk2d_ishock_graph::cyclic_adj_pred(dbsk2d_ishock_edge* e, dbsk2d_ishock_node* v, bool exclude_hidden)
{
  
  if (v->is_a_source()){
    //assume that a source node cannot exist with just one child
    if (e == ((dbsk2d_ishock_node*)v)->cShock())
      return ((dbsk2d_ishock_node*)v)->cShock2();
    else
      return ((dbsk2d_ishock_node*)v)->cShock();
  }
  else if (v->is_an_A3source()){
    //an A3 is always a leaf edge
    return e;
  }
  else {
    //junctions and vcl_sinks can have multiple parents 
    //start looking for the current edge starting from child
    //any one of the edges incident on these nodes can be hidden

    bool bEdgeFound = false;

    dbsk2d_ishock_edge* cSLink = NULL;
    if (v->is_a_junct())
      if ( ((dbsk2d_ishock_node*)v)->cShock())
        if ( ((dbsk2d_ishock_node*)v)->cShock()->isNotHidden() || !exclude_hidden)
          cSLink = ((dbsk2d_ishock_node*)v)->cShock();

    if (e == cSLink)
      bEdgeFound = true;

    //now go through the list of parents
    ishock_edge_list::reverse_iterator curS = v->pShocks().rbegin();
    for(; curS!=v->pShocks().rend(); curS++){
      dbsk2d_ishock_edge* curSLink = (*curS);
      if (curSLink->isHidden() && exclude_hidden) curSLink = NULL;

      if (bEdgeFound && curSLink)
        return curSLink;

      if (e == curSLink)
        bEdgeFound = true;
    }

    if (bEdgeFound){ //this should always be true by now
      if (cSLink)
        return cSLink;
      else {
        //need to go around the parents one more time
        ishock_edge_list::reverse_iterator curS = v->pShocks().rbegin();
        for(; curS!=v->pShocks().rend(); ++curS){
          dbsk2d_ishock_edge* curSLink = (*curS);
          if (curSLink->isNotHidden() || !exclude_hidden) 
            return curSLink;
        }
      }
    }
  }
  
  //it should never get to here
  return NULL;
}


void dbsk2d_ishock_graph::invalid_shocks(
vcl_vector<dbsk2d_ishock_edge*>& shocks)
{


    //first step in validation is to see that the wavefront information
    //exists and is correct
    for ( edge_iterator curE = this->all_edges().begin();
          curE != this->all_edges().end(); curE++ ) 
    {
        dbsk2d_ishock_edge* sedge = (*curE);

        // //make sure that the wavefront structure is still valid
        // if (!sedge->lShock() || !sedge->rShock())
        // {
        //     shocks.push_back(sedge);
        // }

        if ( ob_edges_.count(sedge->id()))
        {
            continue;
        }

        //make sure that the edges terminate at an intersection
        //either at the cell boundary or another edge
        if (!(sedge->cSNode()))
        {
            shocks.push_back(sedge);
        }
    }

}

bool dbsk2d_ishock_graph::valid_shock_graph(bool ignore_ob_shocks)
{

    bool shocks_valid = true;

    //first step in validation is to see that the wavefront information
    //exists and is correct
    for ( edge_iterator curE = this->all_edges().begin();
          curE != this->all_edges().end(); curE++ ) 
    {
        dbsk2d_ishock_edge* sedge = (*curE);

        // //make sure that the wavefront structure is still valid
        // if (!sedge->lShock() || !sedge->rShock())
        // {
        //     shocks_valid = false;
        //     #ifdef DEBUG_SHOCK_VERBOSE
        //     vcl_cout << "S:" << sedge->id() << " wavefronts invalid. \n";
        //     #endif
        //     break;
        // }
        
        if ( !ignore_ob_shocks )
        {
            //make sure that the edges terminate at an intersection
            //either at the cell boundary or another edge
            if (!(sedge->cSNode() || sedge->cell_bnd()))
            {
                shocks_valid = false;
                #ifdef DEBUG_SHOCK_VERBOSE
                vcl_cout << "S:" << sedge->id() << "did not intersect. \n";
                #endif
            }
        }
        else
        {
            if ( !ob_edges_.count(sedge->id()))
            {
                if (!(sedge->cSNode()))
                {
                    shocks_valid = false;
                    #ifdef DEBUG_SHOCK_VERBOSE
                    vcl_cout << "S:" << sedge->id() << "did not intersect. \n";
                    #endif
                }
            }

        }
    }

    //go through the edge_list and purge isolated vertices
    vcl_vector<dbsk2d_ishock_node*> vert_to_del;
    for (vertex_iterator curN = vertices_.begin(); 
         curN != vertices_.end(); curN++){
        dbsk2d_ishock_node* curNode = (*curN);
        if ( curNode->adj_edges().size() == 0 )
        {
            vert_to_del.push_back(curNode);
            break;
        }
       
        
    }

    //is there a better way to do this ???
    for (unsigned int i=0; i<vert_to_del.size(); i++)
    {
        remove_vertex(vert_to_del[i]);
    }
    
    vert_to_del.clear();

    
    #ifdef DEBUG_SHOCK_VERBOSE
    if (shocks_valid)
    {
        vcl_cout << "shocks validated!" << vcl_endl;
    }
    else
    {
        vcl_cout << "shocks invalid!" << vcl_endl;
    }
    #endif

    dbsk2d_assert(shocks_valid);

    if (!shocks_valid){
        vcl_cout << "Shocks computation produced invalid shocks." <<vcl_endl;
    }

    return shocks_valid;
}

void dbsk2d_ishock_graph::update_shocks()
{
  //go through the edge_list and redraw each element that has been flagged
  edge_iterator curE = edges_.begin();
  for (; curE != edges_.end(); curE++){
    dbsk2d_ishock_edge* curShock = (*curE);
    //if (curShock->needs_to_be_redrawn())
      curShock->compute_extrinsic_locus();
  }

  //go through the edge_list and redraw each element that has been flagged
  vertex_iterator curN = vertices_.begin();
  for (; curN != vertices_.end(); curN++){
    dbsk2d_ishock_node* curShock = (*curN);
    //if (curShock->needs_to_be_redrawn())
      curShock->compute_extrinsic_locus();
  }
}

//: Print an ascii summary to the stream
void dbsk2d_ishock_graph::print_summary(vcl_ostream &os)
{
  os << "Shock Graph Summary:" << vcl_endl;
  os << "# of nodes: " << number_of_nodes() << vcl_endl;
  os << "# of edges: " << number_of_edges() << vcl_endl;

  /*os << "E: ";
  
  for (edge_iterator curE = edges_.begin(); curE != edges_.end(); curE++){
    os << (*curE)->id() << " ";
  }
  os << vcl_endl;

  os << "V: ";
  for (vertex_iterator curN = vertices_.begin(); curN != vertices_.end(); curN++)
    os << (*curN)->id() << " ";

  os << vcl_endl;*/
}


void dbsk2d_ishock_graph::ob_shocks()
{

    // clear first
    ob_edges_.clear();

    //first step in validation is to see that the wavefront information
    //exists and is correct
    for ( edge_iterator curE = this->all_edges().begin();
          curE != this->all_edges().end(); curE++ ) 
    {
        dbsk2d_ishock_edge* sedge = (*curE);

        if ( sedge->cell_bnd() )
        {
            ob_edges_[sedge->id()]=sedge;
        }
    }


}

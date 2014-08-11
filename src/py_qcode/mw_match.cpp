/*
This is a code snippet meant for use with scipy.weave.inline.
It relies on the presence of blossom5-v2.04.src/ as a subdir.
If you are using a fresh blossom library, add -fPIC to the Makefile
and then produce the file blossom.o by:
	ld -r *.o -o blossom.o
*/
int edge_idx, vert_idx;
int return_val[num_verts];
PerfectMatching *pm = new PerfectMatching(num_verts, num_edges);

struct PerfectMatching::Options options;
options.verbose = false; //suppress printing from c++

pm->options = options;

for ( edge_idx = 0; edge_idx < num_edges; edge_idx++ )
{
    pm->AddEdge(edges(edge_idx,0), edges(edge_idx,1), edges(edge_idx,2));
}
pm->Solve();
for (vert_idx = 0; vert_idx < num_verts; ++vert_idx)
    {
        int partner = pm->GetMatch(vert_idx);
        partners(vert_idx) = partner;
    }
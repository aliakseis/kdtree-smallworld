// smallworld.cpp : Defines the entry point for the console application.
//

/*
This file is based on a part of ``kdtree'', a library for working with kd-trees:

Copyright (C) 2007-2009 John Tsiombikas <nuclear@siggraph.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <string.h>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>

#include <omp.h>

using std::cerr;
using std::cout;
using std::string;
using std::vector;
using std::ifstream;
using std::exception;
using std::swap;
using std::lower_bound;
using std::nth_element;


/////////////////////////////////////////////////////////////////////////////////

enum { DIM = 2 };

inline double SQ(double x) { return x * x; }

struct kdnode
{
    double pos[DIM];
    int data;
    int dir;
	kdnode *left, *right, *parent;
};

typedef kdnode FriendInfo;


/////////////////////////////////////////////////////////////////////////////////


typedef vector<FriendInfo> FriendInfos;

template<int m_idx>
struct IsInfoLess
{
    bool operator() (const FriendInfo& x, const FriendInfo& y)
    {
        return x.pos[m_idx] < y.pos[m_idx];
    }
    bool operator() (const FriendInfo* x, const FriendInfo* y)
    {
        return x->pos[m_idx] < y->pos[m_idx];
    }
};

template <int dir>
kdnode* insert(vector<FriendInfo*>::iterator begin,
               vector<FriendInfo*>::iterator end,
               kdnode* parent)
{
    if (begin == end)
        return 0;

    int diff = int(end - begin);
    if (1 == diff)
    {
		kdnode* node = *begin;

        node->left = 0;
        node->right = 0;
        node->parent = parent;
        node->dir = dir;

        return node;
    }

    int halfSize = diff / 2;

    vector<FriendInfo*>::iterator middle = begin + halfSize;

	nth_element(begin, middle, end, IsInfoLess<dir>());

	kdnode* node = *middle;

	enum { new_dir = (dir + 1) % DIM };
    node->left = insert<new_dir>(begin, middle, node);
    node->right = insert<new_dir>(++middle, end, node);
    node->parent = parent;
    node->dir = dir;

    return node;
}


struct SearchResult
{
    double dist_sq;
    SearchResult* next;
    int data;
};

bool operator < (const SearchResult& left, const SearchResult& right)
{
    return left.dist_sq < right.dist_sq;
}

bool operator < (const SearchResult& left, double right)
{
    return left.dist_sq < right;
}

bool operator < (double left, const SearchResult& right)
{
    return left < right.dist_sq;
}


class SearchResults
{
public:
    SearchResults()
    : m_pList(0)
    , m_pFree(m_data)
    {
    }

    bool isFull()
    {
        return m_pFree == m_data + 3;
    }
    double dist_sq()
    {
        return m_pList->dist_sq;
    }
    void insert(double dist_sq, kdnode* node)
    {
        if (!isFull())
        {
            SearchResult* newResult = m_pFree++;
            newResult->dist_sq = dist_sq;
            newResult->data = node->data;
            doInsert(&m_pList, newResult);
        }
        else if (m_pList->dist_sq > dist_sq)
        {
            SearchResult* newResult = m_pList;
            m_pList = m_pList->next;
            newResult->dist_sq = dist_sq;
            newResult->data = node->data;
            doInsert(&m_pList, newResult);
        }
    }

    const SearchResult* data()
    {
        return m_pList;
    }

    int size()
    {
        return m_pFree - m_data;
    }


private:
    void doInsert(SearchResult** ppResult, SearchResult* toInsert)
    {
        while ((*ppResult) != 0 && (*ppResult)->dist_sq > toInsert->dist_sq)
            ppResult = &((*ppResult)->next);

        toInsert->next = *ppResult;
        *ppResult = toInsert;
    }

    SearchResult m_data[3];

    SearchResult* m_pList;
    SearchResult* m_pFree;
};


template<int dir>
void kd_nearest_i(kdnode *node, const double *pos,
				  SearchResults& result, double* sq_distances)
{
	kdnode *nearer_subtree, *farther_subtree;

	/* Decide whether to go left or right in the tree */
	const double dist = pos[dir] - node->pos[dir];
	if (dist <= 0) {
		nearer_subtree = node->left;
		farther_subtree = node->right;
	} else {
		nearer_subtree = node->right;
		farther_subtree = node->left;
	}

	enum { new_dir = (dir + 1) % DIM };

	if (nearer_subtree) {
		/* Recurse down into nearer subtree */
		kd_nearest_i<new_dir>(nearer_subtree, pos, result, sq_distances);
	}

    bool do_farther_subtree = false;
    double sq_dist = SQ(dist);
    if (farther_subtree != 0)
    {
		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our minimum distances. */
        if (!result.isFull())
            do_farther_subtree = true;
        else
        {
	        double dist_sq = sq_dist + sq_distances[1 - dir];
            do_farther_subtree = dist_sq < result.dist_sq();
        }
    }

    if (/*node->pos != pos &&*/ (0 == farther_subtree || do_farther_subtree))
    {
	    /* Check the distance of the point at the current node, compare it with our bests so far */
        double dist_sq = sq_dist + SQ(node->pos[1 - dir] - pos[1 - dir]);

        result.insert(dist_sq, node);
    }

	if (do_farther_subtree)
    {
        double save_dist = sq_distances[dir];
        sq_distances[dir] = sq_dist;

		/* Recurse down into farther subtree */
		kd_nearest_i<new_dir>(farther_subtree, pos, result, sq_distances);

        sq_distances[dir] = save_dist;
    }
}


template<int dir>
void kd_nearest_i_nearer_subtree(kdnode *node, const double *pos,
				  SearchResults& result, int* flags)
{
	kdnode *nearer_subtree, *farther_subtree;
    int flag;

	/* Decide whether to go left or right in the tree */
	const double dist = pos[dir] - node->pos[dir];
	if (dist <= 0) {
		nearer_subtree = node->left;
		farther_subtree = node->right;
        flag = 1 << (dir * 2);
	} else {
		nearer_subtree = node->right;
		farther_subtree = node->left;
        flag = 1 << (dir * 2 + 1);
	}

	enum { new_dir = (dir + 1) % DIM };

	if (nearer_subtree) {
		/* Recurse down into nearer subtree */
		kd_nearest_i_nearer_subtree<new_dir>(nearer_subtree, pos, result, flags);
	}

    if (*flags & flag)
        return;

    bool do_farther_subtree = false;
    double sq_dist = SQ(dist);
    if (farther_subtree != 0)
    {
		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our minimum distances. */
        do_farther_subtree = (!result.isFull()
            || sq_dist < result.dist_sq());

        if (!do_farther_subtree)
            *flags |= flag;
    }

    if (node->pos != pos && (0 == farther_subtree || do_farther_subtree))
    {
	    /* Check the distance of the point at the current node, compare it with our bests so far */
        double dist_sq = sq_dist + SQ(node->pos[1 - dir] - pos[1 - dir]);

        result.insert(dist_sq, node);
    }

	if (do_farther_subtree)
    {
        double sq_distances[DIM];// = { 0 };
        sq_distances[1 - dir] = 0;
        sq_distances[dir] = sq_dist;

		/* Recurse down into farther subtree */
		kd_nearest_i<new_dir>(farther_subtree, pos, result, sq_distances);
	}
}


void move_up(kdnode *node, const double *pos,
				  SearchResults& result, int flags)
{
    enum { ALL_SET = (1 << (DIM * 2)) - 1 };

    while (flags != ALL_SET)
    {
        kdnode* childNode = node;
        node = node->parent;
        if (0 == node)
        {
            break;
        }

        kdnode *farther_subtree;
        int flag;

        const int dir = node->dir;

        // Decide whether to go left or right in the tree
        //if (dist <= 0) {
        if (childNode == node->left)
        {
            flag = 1 << (dir * 2);
            if (flags & flag)
                continue;
            farther_subtree = node->right;
        } else {
            flag = 1 << (dir * 2 + 1);
            if (flags & flag)
                continue;
            farther_subtree = node->left;
        }

        const double dist = pos[dir] - node->pos[dir];


        bool do_farther_subtree = false;
        double sq_dist = SQ(dist);
        if (farther_subtree != 0)
        {
            // Check if we have to recurse down by calculating the closest
            // point of the hyperrect and see if it's closer than our minimum distances.
            do_farther_subtree = (!result.isFull()
                || sq_dist < result.dist_sq());

            if (!do_farther_subtree)
                flags |= flag;
        }

        if (0 == farther_subtree || do_farther_subtree)
        {
            // Check the distance of the point at the current node, compare it with our bests so far
            double dist_sq = sq_dist + SQ(node->pos[1 - dir] - pos[1 - dir]);

            result.insert(dist_sq, node);
        }

        if (do_farther_subtree)
        {
            double sq_distances[DIM];// = { 0 };
            // Recurse down into farther subtree
            switch (dir)
            {
            case 0:
                sq_distances[1] = 0;
                sq_distances[0] = sq_dist;//SQ(dist);
                kd_nearest_i<1>(farther_subtree, pos, result, sq_distances);
                break;
            case 1:
                sq_distances[0] = 0;
                sq_distances[1] = sq_dist;//SQ(dist);
                kd_nearest_i<0>(farther_subtree, pos, result, sq_distances);
                break;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////

typedef int OneResult[4];

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cerr << "Usage: smallworld input_file\n";
		return 1;
	}

	try
	{
	    /*
		string path(argv[1]);
		if (string::npos == path.find_first_of("\\/"))
		{
			string dir(argv[0]);
			string::size_type pos = dir.find_last_of("\\/");
			if (string::npos != pos)
			{
				path = dir.substr(0, pos + 1) + path;
			}
		}

		ifstream inFile(path.c_str());
		//*/

        FILE* f = fopen(argv[1], "r");
        FriendInfos infos;

        FriendInfo info;
		while (3 == fscanf(f, "%d %le %le", &info.data, &info.pos[0], &info.pos[1]))
        {
            infos.push_back(info);
        }

        fclose(f);

        if (infos.size() < 2)
            return 0;

        vector<FriendInfo*> infoPtrs;
        infoPtrs.reserve(infos.size());

        for (FriendInfos::iterator it = infos.begin(); it != infos.end(); ++it)
        {
            infoPtrs.push_back(&*it);
        }

        kdnode* root = insert<0>(infoPtrs.begin(), infoPtrs.end(), 0);


        OneResult* results = new OneResult[infos.size()];

        const int numInfos = infos.size();
        FriendInfo* const arrInfos = &infos[0];
        FriendInfo** const arrInfoPtrs = &infoPtrs[0];

        int i;
#pragma omp parallel for schedule (dynamic, 10)
        for (i = 0; i < numInfos; ++i)
        {
            SearchResults result;

            int flags = 0;

            /* Search for the nearest neighbours recursively */
            //kdnode* node = arrInfos + i;
            kdnode* node = arrInfoPtrs[i];

            switch (node->dir)
            {
            case 0:
                kd_nearest_i_nearer_subtree<0>(node, node->pos, result, &flags);
                break;
            case 1:
                kd_nearest_i_nearer_subtree<1>(node, node->pos, result, &flags);
                break;
            }

            move_up(node, node->pos, result, flags);

            //OneResult& out = results[i];
            OneResult& out = results[node - arrInfos];

            out[0] = node->data;

            const SearchResult* pResult =  result.data();
            for (int j = (int) result.size(); --j >= 0;)
            {
                out[j + 1] = pResult->data;
                pResult = pResult->next;
            }
        }

        for (i = 0; i < (int) infos.size(); ++i)
        {
            const OneResult& result = results[i];

            switch (infos.size())
            {
            case 3:
                printf("%d %d,%d\n", result[0], result[1], result[2]);
                break;
            case 2:
                printf("%d %d\n", result[0], result[1]);
                break;
            default:
                printf("%d %d,%d,%d\n", result[0], result[1], result[2], result[3]);
                break;
            }
        }

        delete[] results;
	}
	catch (exception& e)
	{
		cerr << e.what();
		return 1;
	}

	return 0;
}


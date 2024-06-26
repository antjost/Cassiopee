#include <algorithm>
#include <cassert>

#include "hedge.h"
#include "vertex.h"
#include "primitives.h"
#include "dcel.h"

Hedge::Hedge(Vertex *v)
: orig(v), twin(NULL), prev(NULL), next(NULL), left(NULL), color(Dcel::NO_IDEA),
  cycle(NULL)
{}

static
E_Int _partition(std::vector<Hedge *> &H, E_Int low, E_Int high)
{
    Hedge *pivot = H[high];
    E_Int i = low-1;

    for (E_Int j = low; j < high; j++) {
        if (Hedge::cmp_cwise(H[j], pivot) <= 0) {
            i++;
            std::swap(H[i], H[j]);
        }
    }

    i++;
    std::swap(H[i], H[high]);
    return i;
}

void Hedge::sort_cwise(std::vector<Hedge *> &H, E_Int low, E_Int high)
{
    if (low >= high)
        return;

    E_Int p = _partition(H, low, high);

    sort_cwise(H, low, p - 1);
    sort_cwise(H, p + 1, high);
}

void Hedge::sort_ccwise(std::vector<Hedge *> &H, E_Int low, E_Int high)
{
    sort_cwise(H, low, high);
    std::reverse(H.begin(), H.end());
}

E_Int Hedge::cmp_cwise(const Hedge *h, const Hedge *w)
{
    assert(h->orig == w->orig);
    Vertex *c = h->orig;
    Vertex *a = h->twin->orig;
    Vertex *b = w->twin->orig;

    E_Float ax = a->x;
    E_Float ay = a->y;
    E_Float bx = b->x;
    E_Float by = b->y;
    E_Float cx = c->x;
    E_Float cy = c->y;

    long double acx = (long double)ax - (long double)cx;
    long double acy = (long double)ay - (long double)cy;
    long double bcx = (long double)bx - (long double)cx;
    long double bcy = (long double)by - (long double)cy;

    E_Int sign_acx = Sign(acx);
    E_Int sign_acy = Sign(acy);
    E_Int sign_bcx = Sign(bcx);
    E_Int sign_bcy = Sign(bcy);

    if (sign_acx >= 0 && sign_bcx < 0)
        return -1;
    if (sign_acx < 0 && sign_bcx >= 0)
        return 1;
    if (sign_acx == 0 && sign_bcx == 0) {
        if (sign_acy >= 0 || sign_bcy >= 0) {
            long double diff = (long double)ay - (long double)by;
            if (Sign(diff) > 0) return -1;
            else return 1;
        }

        long double diff = (long double)by - (long double)ay;
        if (Sign(diff) > 0) return -1;
        else return 1;
    }
    
    E_Float det = acx * bcy - bcx * acy;
    //E_Float det = DifferenceOfProducts(acx, bcy, bcx, acy);
    E_Int cmp = Sign(det);

    if (cmp < 0)
        return -1;
    else if (cmp > 0)
        return 1;

    // Overlapping segments
    
    assert(h->color != w->color);

    // If right half, red before black
    // Otherwise, black before red

    cmp = Sign(h->color - w->color);

    if (sign_acx >= 0) {
        assert(sign_bcx >= 0);
        return cmp;
    } else {
        return -cmp;
    }
}
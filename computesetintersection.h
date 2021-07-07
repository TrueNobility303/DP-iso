#ifndef SUBGRAPHMATCHING_COMPUTE_SET_INTERSECTION_H
#define SUBGRAPHMATCHING_COMPUTE_SET_INTERSECTION_H

#include "types.h"
#include "config.h"
#include <immintrin.h>
#include <x86intrin.h>
#include <cstdint>

/*
 * Because the set intersection is designed for computing common neighbors, the target is uieger.
 */

class ComputeSetIntersection {
public:

    static void ComputeCandidates(const VertexID* larray, ui l_count, const VertexID* rarray,
                                  ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCandidates(const VertexID* larray, ui l_count, const VertexID* rarray,
                                  ui r_count, ui &cn_count);

    static void ComputeCNNaiveStdMerge(const VertexID* larray, ui l_count, const VertexID* rarray,
                                       ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCNNaiveStdMerge(const VertexID* larray, ui l_count, const VertexID* rarray,
                                       ui r_count, ui &cn_count);

    static void ComputeCNGalloping(const VertexID * larray, ui l_count, const VertexID * rarray,
                                   ui r_count, VertexID * cn, ui& cn_count);
    static void ComputeCNGalloping(const VertexID * larray, ui l_count, const VertexID * rarray,
                                   ui r_count, ui& cn_count);
    static const ui GallopingSearch(const VertexID *src, ui begin, ui end, ui target);
    static const ui BinarySearch(const VertexID *src, ui begin, ui end, ui target);

};


void ComputeSetIntersection::ComputeCandidates(const VertexID* larray, const ui l_count,
                                               const VertexID* rarray, const ui r_count,
                                               VertexID* cn, ui &cn_count) {
    return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn, cn_count);
}

void ComputeSetIntersection::ComputeCandidates(const VertexID* larray, const ui l_count,
                                               const VertexID* rarray, const ui r_count,
                                               ui &cn_count) {
    return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn_count);
}

//naive方法集合求交，遍历两个数组记录元素求交集
void ComputeSetIntersection::ComputeCNNaiveStdMerge(const VertexID* larray, const ui l_count,
                                                    const VertexID* rarray, const ui r_count,
                                                    VertexID* cn, ui &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    ui lc = l_count;
    ui rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        ui tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    ui li = 0;
    ui ri = 0;

    while (true) {
        if (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }
        else if (larray[li] > rarray[ri]) {
            ri += 1;
            if (ri >= rc) {
                return;
            }
        }
        else {
            cn[cn_count++] = larray[li];

            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNNaiveStdMerge(const VertexID* larray, const ui l_count,
                                                    const VertexID* rarray, const ui r_count,
                                                    ui &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    ui lc = l_count;
    ui rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        ui tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    ui li = 0;
    ui ri = 0;

    while (true) {
        if (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }
        else if (larray[li] > rarray[ri]) {
            ri += 1;
            if (ri >= rc) {
                return;
            }
        }
        else {
            cn_count += 1;
            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNGalloping(const VertexID* larray, const ui l_count,
                                                const VertexID* rarray, const ui r_count,
                                                VertexID* cn, ui &cn_count) {
    ui lc = l_count;
    ui rc = r_count;
    cn_count = 0;
    if (lc == 0 || rc == 0)
        return;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        ui tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    ui li = 0;
    ui ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = GallopingSearch(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn[cn_count++] = larray[li];

            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNGalloping(const VertexID* larray, const ui l_count,
                                                const VertexID* rarray, const ui r_count,
                                                ui &cn_count) {
    ui lc = l_count;
    ui rc = r_count;
    cn_count = 0;
    if (lc == 0 || rc == 0)
        return;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        ui tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    ui li = 0;
    ui ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = GallopingSearch(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn_count += 1;

            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

const ui ComputeSetIntersection::GallopingSearch(const VertexID *src, const ui begin, const ui end,
                                            const ui target) {
    if (src[end - 1] < target) {
        return end;
    }
    // galloping
    if (src[begin] >= target) {
        return begin;
    }
    if (src[begin + 1] >= target) {
        return begin + 1;
    }
    if (src[begin + 2] >= target) {
        return begin + 2;
    }

    ui jump_idx = 4;
    ui offset_beg = begin;
    while (true) {
        ui peek_idx = offset_beg + jump_idx;
        if (peek_idx >= end) {
            return BinarySearch(src, (jump_idx >> 1) + offset_beg + 1, end, target);
        }
        if (src[peek_idx] < target) {
            jump_idx <<= 1;
        } else {
            return src[peek_idx] == target ? peek_idx :
                   BinarySearch(src, (jump_idx >> 1) + offset_beg + 1, peek_idx + 1, target);
        }
    }
}

const ui ComputeSetIntersection::BinarySearch(const VertexID *src, const ui begin, const ui end, const ui target) {
    int offset_begin = begin;
    int offset_end = end;
    while (offset_end - offset_begin >= 16) {
        auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_begin) + offset_end) / 2);
        _mm_prefetch((char *) &src[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &src[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
        if (src[mid] == target) {
            return mid;
        } else if (src[mid] < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (src[offset] >= target) {
            return (ui)offset;
        }
    }

    return (ui)offset_end;
}

#endif //FSE_COMPUTESETINTERSECTION_H
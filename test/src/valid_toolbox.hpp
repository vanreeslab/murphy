#include <list>

#include "core/types.hpp"
#include "tools/patch.hpp"

static void TreeCaseIdToPatch(const level_t start_level, const int case_id, std::list<Patch>* patch_list) {
    for (int itree = 0; itree < 8; ++itree) {
        real_t origin[3] = {1.0 * (itree % 2),
                            1.0 * ((itree % 4) / 2),
                            // 1.0 * (itree % 4) / 2,
                            1.0 * (itree / 4)};
        real_t length[3] = {1.0, 1.0, 1.0};

        level_t level = (case_id >> itree) % 2 + start_level;
        patch_list->push_back(Patch(origin, length, level));

        m_log("tree %d @ %f %f %f has level %d", itree, origin[0], origin[1], origin[2], level);
    }
}
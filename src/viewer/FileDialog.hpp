//
// Created by Lei on 4/4/2024.
//

#ifndef PCO_FILEDIALOG_HPP
#define PCO_FILEDIALOG_HPP

#include "Config.hpp"

#include <cstdio>
#include <string>

#ifdef _WIN32

#include <windows.h>
#include <Commdlg.h>

#endif

NAMESPACE_BEGIN(PCO)
    namespace viewer {

        std::string fileDialogOpen();

        std::string fileDialogSave();

    } // namespace viewer
NAMESPACE_END(PCO)

#endif //PCO_FILEDIALOG_HPP

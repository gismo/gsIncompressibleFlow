/** @file gsFlowLogger.h
 
    @brief Logger class for incompressible flow solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gismo.h>

using namespace gismo;

class gsFlowLogger
{

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsFlowLogger> Ptr;
    typedef memory::unique_ptr<gsFlowLogger> uPtr;


public: // *** Enums ***

    enum class mode { terminal, file, all, quiet };


protected: // *** Class members ***

    int m_rank;
    mode m_mode;
    std::unique_ptr<std::ofstream> m_file;


public: // *** Constructor/destructor ***

    gsFlowLogger(mode mode, const std::string& filename = "output.log", int rank = 0)
        : m_mode(mode), m_rank(rank)
    {
        if (m_mode == mode::file || m_mode == mode::all)
        {
            m_file = std::make_unique<std::ofstream>(filename);

            if (!m_file->is_open())
                GISMO_ERROR("Cannot open log file: " + filename);
        }
    }

public: // *** Static functions ***

    static mode parseOutputMode(const std::string& s)
    {
        if (s == "terminal") return mode::terminal;
        if (s == "file")     return mode::file;
        if (s == "all")      return mode::all;
        if (s == "quiet")    return mode::quiet;

        gsWarn << "Unknown output mode: '" << s << "', using 'terminal' instead.\n";
        return mode::terminal;
    }


public: // *** Member functions ***

    template <typename T>
    gsFlowLogger& log(const T& msg, bool fileOnly = false)
    {
        if (m_rank != 0) return *this; // only rank 0 logs

        if (!fileOnly && (m_mode == mode::terminal || m_mode == mode::all))
            gsInfo << msg;

        if ((m_mode == mode::file || m_mode == mode::all) && m_file)
            (*m_file) << msg;

        return *this;
    }

    template <typename T>
    gsFlowLogger& operator<<(const T& msg) { return log(msg); }


public: // *** Getters/setters ***

    mode getMode() const { return m_mode; }
    void setMode(mode mode) { m_mode = mode; }
    
};
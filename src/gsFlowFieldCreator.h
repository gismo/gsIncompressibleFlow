/** @file gsFlowFieldCreator.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once

namespace gismo
{

/**
 * @brief Function returning the (unit) outer normal vector scaled by values of a given scalar field.
 * 
 * The function is defined on a particular side of a particular patch.
 * Useful for definition of boundary conditions for flow problems, e.g. in the form pressure * normal.
 * 
 * @tparam T real number type
 */
template<class T>
class gsScaledOuterNormalField : public gsFunction<T> 
{

public: // *** Smart pointers ***

    typedef memory::shared_ptr< gsScaledOuterNormalField > Ptr;
    typedef memory::unique_ptr< gsScaledOuterNormalField > uPtr;


protected: // *** Class members ***

    index_t m_patch;
    boxSide m_side;
    const gsField<T>& m_scalarField;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] patch        patch index
    /// @param[in] side         patch side
    /// @param[in] scalarField  the scalar field for scaling
    gsScaledOuterNormalField( index_t patch, boxSide side, const gsField<T>& scalarField):
    m_patch(patch), m_side(side), m_scalarField(scalarField)
    { }

    GISMO_CLONE_FUNCTION(gsScaledOuterNormalField)


public: // *** Member functions ***

    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(this->targetDim(), u.cols());

        gsMapData<T> mapData(NEED_OUTER_NORMAL);
        mapData.side = m_side;
        mapData.points = u;
        m_scalarField.patches().patch(m_patch).computeMap(mapData);

        gsMatrix<T> scalarVals = getScalarValues(u);

        gsVector<> normal;
        for (index_t k = 0; k < u.cols(); k++)
        {
            outerNormal(mapData, k, m_side, normal);
            result.col(k) = scalarVals(0, k) * normal;
        }
    }

    short_t domainDim() const { return m_scalarField.patches().patch(m_patch).parDim(); }
    short_t targetDim() const { return m_scalarField.patches().patch(m_patch).geoDim(); }

    virtual std::ostream &print(std::ostream &os) const
    { os << "gsScaledOuterNormalField"; return os; };


protected: // *** Member functions ***

    virtual gsMatrix<T> getScalarValues(const gsMatrix<T>& u) const
    { return m_scalarField.function(m_patch).eval(u); }


}; // gsScaledOuterNormalField

// ----------------------------------------------------------------------------

/**
 * @brief Function returning the (unit) outer normal vector scaled by the difference of a target value and values of a given scalar field.
 * 
 * The function is defined on a particular side of a particular patch.
 * Useful for definition of boundary conditions for flow problems, e.g. in the form (targetPressure - pressure) * normal.
 * 
 * @tparam T real number type
 */
template<class T>
class gsDiffScaledOuterNormalField : public gsScaledOuterNormalField<T> 
{

public:
    typedef gsScaledOuterNormalField<T> Base;

public: // *** Smart pointers ***

    typedef memory::shared_ptr< gsDiffScaledOuterNormalField > Ptr;
    typedef memory::unique_ptr< gsDiffScaledOuterNormalField > uPtr;


protected: // *** Class members ***

    T m_targetValue;
    

public: // *** Constructor/destructor ***

    gsDiffScaledOuterNormalField( index_t patch, boxSide side, const gsField<T>& scalarField, T targetValue):
    Base(patch, side, scalarField), m_targetValue(targetValue)
    { }

    GISMO_CLONE_FUNCTION(gsDiffScaledOuterNormalField)


public: // *** Member functions ***

    virtual std::ostream &print(std::ostream &os) const
    { os << "gsDiffScaledOuterNormalField"; return os; };


protected: // *** Member functions ***

    virtual gsMatrix<T> getScalarValues(const gsMatrix<T>& u) const
    { 
        gsMatrix<T> fieldVals = Base::getScalarValues(u);

        gsMatrix<T> targetVals(fieldVals.rows(), fieldVals.cols());
        targetVals.setOnes();
        targetVals *= m_targetValue;

        return (targetVals - fieldVals); 
    }


}; // gsDiffScaledOuterNormalField

} // namespace gismo
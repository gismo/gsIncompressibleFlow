/** @file uwbObjectiveFunctionEvaluator.h

    Author(s): E. Turnerova

*/

#pragma once

namespace gismo
{

template <class T>
class uwbObjectiveFunctionEvaluator
{
public:
    uwbObjectiveFunctionEvaluator()
    {
        initMembers();
    }

    uwbObjectiveFunctionEvaluator(const gsField<T>& uSolField,
                                  const gsField<T>& pSolField,
                                  const std::vector< gsMultiBasis<T> >& bases,
                                  const gsMultiPatch<T>& patches) :
        m_uSolField(uSolField), m_pSolField(pSolField), m_bases(bases), m_patches(patches)
    {
        initMembers();
    }

    uwbObjectiveFunctionEvaluator(const gsField<T>& loftedPressure) : m_loftedPressure(loftedPressure){ initMembers(); }

    virtual ~uwbObjectiveFunctionEvaluator()
    { }

    void initMembers()
    {
        m_membersReady = false;
        m_objectiveEvaluated = false;
        m_objectiveAbsolute = true;

        m_weightLift = 0.0;
        m_lift = 0.0;
        m_lift0 = 1.0;
        m_radius = 0.0;
        m_weightOutflowVelocity = 0.0;
        m_outflowVelocity = 0.0;
        m_outflowVelocity0 = 1.0;
        m_averageVel = 0.0;
        m_uTangentialTarget = 0.0;
        m_velAvg = false;
        m_weightPressureProfile = 0.0;
        m_pressureProfile = 0.0;
        m_pressureProfile0 = 1.0;
        m_pTargetProfile.setZero(2);
        m_weightHead = 0.0;
        m_head = 0.0;
        m_targetHead = 0.0;
        m_weightEfficiency = 0.0;
        m_efficiency = 0.0;
        m_efficiency0 = 1.0;
        m_weightCavitation = 0.0;
        m_cavitation = 0.0;

        m_tarDim = 0.0;

        m_angularVel = 0.;
        m_flowRate = 0.;
        m_gravity = 0.;
        m_numBlades = 0.;
        m_density = 0.;
    }

    void initialize(gsVector<int> profilePatchNumbers, std::vector<boxSide> profileSides, gsVector<T> normalStream, T radius,
                    gsVector<int> outPatchNumbers, std::vector<boxSide> outSides, T uTangentialTarget, bool velAvg = false)
    {
        m_normalStream = normalStream;
        m_profileSides = profileSides;
        m_profilePatches = profilePatchNumbers;
        m_radius = radius;
        m_outSides = outSides;
        m_outPatches = outPatchNumbers;
        m_uTangentialTarget = uTangentialTarget;
        m_velAvg = velAvg;

        m_tarDim = m_patches.dim();

        m_membersReady = true;
    }

    void initialize(gsVector<int> bladePatchNumbers, std::vector<boxSide> bladeSides,
                    gsVector<int> outPatchNumbers, std::vector<boxSide> outSides, T uTangentialTarget,
                    gsVector<int> inPatchNumbers, std::vector<boxSide> inSides,
                    T head)
    {
        m_bladeSides = bladeSides;
        m_bladePatches = bladePatchNumbers;
        m_outSides = outSides;
        m_outPatches = outPatchNumbers;
        m_uTangentialTarget = uTangentialTarget;
        m_inSides = inSides;
        m_inPatches = inPatchNumbers;
        m_head = head;

        m_tarDim = m_patches.dim();

        m_membersReady = true;
    }

    void initialize(T angularVel, T flowRate, T gravity, int numBlades, short_t tarDim, T density)
    {
        m_angularVel = angularVel;
        m_flowRate = flowRate;
        m_gravity = gravity;
        m_numBlades = numBlades;
        m_density = density;

        m_tarDim = tarDim;

        m_membersReady = true;
    }

    void updateSolution(const gsField<T>& uSolField, const gsField<T>& pSolField)
    {
        m_uSolField = uSolField;
        m_pSolField = pSolField;
    }

    void updateSolution(const gsField<T>& loftedPressure) { m_loftedPressure = loftedPressure; }
    //void updateBasis(const gsBasis<T>& basis) { m_basis = basis; }

public:

    T evaluate()
    {
        GISMO_ENSURE(checkSumOfWeights() == true && m_membersReady == true, "Sum of weights does not match or not all members are initialized.");

        switch (m_tarDim)
        {
        case 2:
            m_lift = 0.0;
            m_outflowVelocity = 0.0;
            m_pressureProfile = 0.0;
            for (int k = 0; k < m_profilePatches.rows(); k++)
            {
                computeLift(k);
                computePressureProfilePart(k);
            }

            for (int k = 0; k < m_outPatches.rows(); k++)
            {
                if (m_velAvg)
                    computeOutflowVelocityAvg(k, 1); //second input parameter represents component of the velocity which is averaged - in 2D it is 'y' component as tangential velocity
                else
                    computeOutflowVelocity(k);
            }

            m_objectiveEvaluated = true;
            return (m_weightLift * getLift() + m_weightOutflowVelocity * getOutflowVelocity() + m_weightPressureProfile * getPressureProfile());
            break;
        case 3:
            m_head = 0.0;
            m_efficiency = 0.0;
            m_outflowVelocity = 0.0;
            m_cavitation = 0.0;
            computeHead();
            computeEfficiency();
            for (int k = 0; k < m_outPatches.rows(); k++)
                computeOutflowVelocity(k);
            computeCavitation();

            m_objectiveEvaluated = true;
            return (m_weightHead * m_head + m_weightEfficiency * m_efficiency + m_weightOutflowVelocity * m_outflowVelocity + m_weightCavitation * m_cavitation);
            break;
        default:
            GISMO_ERROR("Optimization implemented only for 2D and 3D.");
            return 0;
            break;
        }
    }

    void computeLift(int index)
        //computes integral (p*n_p*n_s) over patch side, where p is pressure, n_p is unit outer normal of the patch side, n_s is the normal vector to the stream
    {
        const gsBasis<T> & basis = m_bases.at(1).basis(m_profilePatches[index]);

        gsVector<int> numQuadNodes(m_tarDim);
        const int dir = m_profileSides[index].direction();
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = (2 * basis.degree(i) + 1);
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here

        gsMatrix<T> quNodes; // Mapped nodes
        gsVector<T> quWeights; // Mapped weights
        unsigned evFlags = NEED_VALUE | NEED_OUTER_NORMAL;   // Evaluation flags for the Geometry map

        // Initialize geometry evaluator
         typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_patches.patch(m_profilePatches[index])));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(m_profileSides[index]);
        for (; domIt->good(); domIt->next())
        {
            // Compute the quadrature rule
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval->evaluateAt(quNodes);

            // Evaluate solution on element nodes
            gsMatrix<T> solPVals = m_pSolField.value(quNodes, m_profilePatches[index]);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Compute the outer normal vector
                gsVector<T> normal;
                geoEval->outerNormal(k, m_profileSides[index], normal);
                m_lift += quWeights[k] * solPVals(k) * normal.dot(m_normalStream.transpose());

            }
        }
    }

    void computeInitialLift()
    {
        for (int k = 0; k < m_profilePatches.rows(); k++)
            computeLift(k);

        m_lift0 = getLift();
        m_lift = 0.0;
    }

    //ready for 2D case, where the axis of rotation equals x-axis
    //WILL BE GENERALIZED!!!
    void computeOutflowVelocity(int index)
    {
        switch (m_tarDim)
        {
        case 2:
            computeOutflowVelocity2D(index);
            break;
        case 3:
            computeOutflowVelocity3D(index);
            break;
        default:
            GISMO_ERROR("Optimization implemented only for 2D and 3D.");
            break;
        }
    }

    void computeOutflowVelocity2D(int index)
    {
        const gsBasis<T> & basis = m_bases.at(0).basis(m_outPatches[index]);

        gsVector<int> numQuadNodes(m_tarDim);
        const int dir = m_outSides[index].direction();
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = (2 * basis.degree(i) + 1);
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here

        gsMatrix<T> quNodes; // Mapped nodes
        gsVector<T> quWeights; // Mapped weights
        unsigned evFlags = NEED_VALUE | NEED_OUTER_NORMAL;   // Evaluation flags for the Geometry map

        // Initialize geometry evaluator
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_patches.patch(m_outPatches[index])));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(m_outSides[index]);
        for (; domIt->good(); domIt->next())
        {
            // Compute the quadrature rule on patch1
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval->evaluateAt(quNodes);

            // Evaluate solution on element nodes
            gsMatrix<T> solUVals = m_uSolField.value(quNodes, m_outPatches[index]);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                gsVector<T> normal;
                geoEval->outerNormal(k, m_outSides[index], normal);

                // weight * abs(det J), where J is geometry Jacobian
                // geoEval->measure(k) and geoEval->jacDet(k) not working as geoEval is initialized at the patch, not on boundary
                const T weight = quWeights[k] * normal.norm();

                //tangential velocity equals u_y in 2D profile example
                m_outflowVelocity += weight * math::pow(solUVals(1, k) - m_uTangentialTarget, 2);
            }
        }
    }

    void computeOutflowVelocity3D(int index)
    {
        gsWarn << "computeOutflowVelocity() for 3D not implemented yet.\n";
    }

    T evaluateOutflowVelocityAvg(int index, int component)
    {
        computeOutflowVelocityAvg(index, component);
        return m_averageVel;
    }

    void computeOutflowVelocityAvg(int index, int component)
    {
        GISMO_ASSERT(component >= 0 && component < m_tarDim, "Invalid velocity component index.");

        const gsBasis<T> & basis = m_bases.at(0).basis(m_outPatches[index]);

        gsVector<int> numQuadNodes(m_tarDim);
        const int dir = m_outSides[index].direction();
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = (2 * basis.degree(i) + 1);
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here

        gsMatrix<T> quNodes; // Mapped nodes
        gsVector<T> quWeights; // Mapped weights
        unsigned evFlags = NEED_VALUE | NEED_OUTER_NORMAL;   // Evaluation flags for the Geometry map

                                                             // Initialize geometry evaluator
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_patches.patch(m_outPatches[index])));

        T area = 0.0;
        T velIntegral = 0.0;

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(m_outSides[index]);

        for (; domIt->good(); domIt->next())
        {
            // Compute the quadrature rule on patch1
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval->evaluateAt(quNodes);

            // Evaluate solution on element nodes
            gsMatrix<T> solUVals = m_uSolField.value(quNodes, m_outPatches[index]);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                gsVector<T> normal;
                geoEval->outerNormal(k, m_outSides[index], normal);

                const T weight = quWeights[k] * normal.norm();

                //tangential velocity equals u_y in 2D profile example
                area += weight;
                velIntegral += weight * solUVals(component, k);
            }
        }
        m_averageVel = (velIntegral / area);
        m_outflowVelocity = math::abs(m_averageVel - m_uTangentialTarget);//math::pow(m_averageVel - m_uTangentialTarget, 2);
    }

    void computePressureProfilePart(int index)
    {
        const gsBasis<T> & basis = m_bases.at(1).basis(m_profilePatches[index]);

        gsVector<int> numQuadNodes(m_tarDim);
        const int dir = m_profileSides[index].direction();
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = (2 * basis.degree(i) + 1);
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here

        gsMatrix<T> quNodes; // Mapped nodes
        gsVector<T> quWeights; // Mapped weights
        unsigned evFlags = NEED_VALUE | NEED_MEASURE;   // Evaluation flags for the Geometry map

        // Initialize geometry evaluator
       typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_patches.patch(m_profilePatches[index])));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(m_profileSides[index]);
        for (; domIt->good(); domIt->next())
        {
            // Compute the quadrature rule
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval->evaluateAt(quNodes);

            // Evaluate solution on element nodes
            gsMatrix<T> solPVals = m_pSolField.value(quNodes, m_profilePatches[index]);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                gsVector<T> normal;
                geoEval->outerNormal(k, m_profileSides[index], normal);

                // weight * abs(det J), where J is geometry Jacobian
                // geoEval->measure(k) and geoEval->jacDet(k) not working as geoEval is initialized at the patch, not on boundary
                const T weight = quWeights[k] * normal.norm();
                m_pressureProfile += weight * math::pow(solPVals(k) - m_pTargetProfile[index], 2);
            }
        }
    }

    T evaluateEfficiencyFromLoftedPressure(const gsBasis<T>& basis, int dir = 1)
    {
        GISMO_ENSURE(m_weightEfficiency > 0 && m_membersReady == true && m_head > 0, "Not all members are initialized, weight of efficiency or head is not set.");

        m_efficiency = 0.;
        T torqueMoment = 0.;

        const gsGeometry<T>& geo = m_loftedPressure.geometry();

        gsVector<index_t> numQuadNodes(m_tarDim);
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = 2 * basis.maxDegree() + 1;

        gsGaussRule<T> QuRule(numQuadNodes);

        gsMatrix<T> quNodes;
        gsVector<T> quWeights;

        unsigned evFlags = NEED_VALUE | NEED_OUTER_NORMAL | NEED_MEASURE;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, geo));//(patches.patch(0).evaluator(evFlags));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            geoEval->evaluateAt(quNodes);
            gsMatrix<T> solPVals = m_loftedPressure.value(quNodes);
            gsMatrix<T> physNodes, guideNodes, guideVec;
            geo.eval_into(quNodes, physNodes); // physical coordinates
            guideNodes.setZero(physNodes.rows(), physNodes.cols());
            guideNodes.row(0) = physNodes.row(0);
            guideVec = physNodes - guideNodes;

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Compute the outer normal vector
                gsVector<T> normal;
                geoEval->normal(k, normal);
                normal /= dir * normal.norm();

                const T weight = quWeights(k) * geoEval->measure(k);
                torqueMoment += weight * solPVals(k) * m_density * (normal(2) * guideVec(1,k) - normal(1) * guideVec(2,k));
            }
        }
        torqueMoment *= m_numBlades;
        gsInfo << "torqueMoment = " << torqueMoment << "\n";
        m_efficiency = torqueMoment * m_angularVel / (m_flowRate * m_density * m_gravity * m_head);

        if (isObjectiveAbsolute())
            return m_weightEfficiency * m_efficiency;
        else
            return m_weightEfficiency * m_efficiency0 / m_efficiency;
    }

    T evaluateHeadFromLoftedPressure(const gsBasis<T>& basis, const gsField<T> loftedPressure, T velMeridial, T gravity)
    {
        T totalPressure = 0.;
        T area = 0.;

        const gsGeometry<T>& geo = loftedPressure.geometry();

        gsVector<index_t> numQuadNodes(m_tarDim);
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = 2 * basis.maxDegree() + 1;

        gsGaussRule<T> QuRule(numQuadNodes);

        gsMatrix<T> quNodes;
        gsVector<T> quWeights;
        unsigned evFlags = NEED_VALUE | NEED_MEASURE;
       typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, geo));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            geoEval->evaluateAt(quNodes);
            gsMatrix<T> solPVals = loftedPressure.value(quNodes);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                const T weight = quWeights(k) * geoEval->measure(k);
                area += weight;
                totalPressure += weight * (solPVals(k) + 0.5 * math::pow(velMeridial, 2));
            }
        }
        totalPressure /= (area * gravity);

        gsInfo << "area = " << area << "\n";

    return totalPressure;
    }

    T evaluateHeadFromLoftedPressure(const gsField<T> loftedPressure, const gsBasis<T>& basis, const gsField<T> loftedVelocity, T gravity)
    {
        T totalPressure = 0.;
        T area = 0.;

        const gsGeometry<T>& geo = loftedVelocity.geometry();

        gsVector<index_t> numQuadNodes(m_tarDim);
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = 2 * basis.maxDegree() + 1;

        gsGaussRule<T> QuRule(numQuadNodes);

        gsMatrix<T> quNodes;
        gsVector<T> quWeights;
        unsigned evFlags = NEED_VALUE | NEED_MEASURE | NEED_NORMAL;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, geo));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            geoEval->evaluateAt(quNodes);
            gsMatrix<T> solPVals = loftedPressure.value(quNodes);
            gsMatrix<T> solUVals = loftedVelocity.value(quNodes);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                gsVector<T> normal;
                geoEval->normal(k, normal);
                normal /= normal.norm();

                const T weight = quWeights(k) * geoEval->measure(k);
                area += weight;
                //totalPressure += weight * (solPVals(k) + 0.5 * math::pow(solUVals(0, k), 2));
                totalPressure += weight * (solPVals(k) + 0.5 * math::pow(normal.dot(solUVals.col(k)), 2));
            }
        }
        totalPressure /= (area * gravity);

        gsInfo << "area = " << area << "\n";

        return totalPressure;
    }

    T evaluateFlowRateFromLoftedVelocity(const gsBasis<T>& basis, const gsField<T> loftedVelocity) const
    {
        T flowRate = 0.;

        const gsGeometry<T>& geo = loftedVelocity.geometry();

        gsVector<index_t> numQuadNodes(m_tarDim);
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = 2 * basis.maxDegree() + 1;

        gsGaussRule<T> QuRule(numQuadNodes);

        gsMatrix<T> quNodes;
        gsVector<T> quWeights;
        unsigned evFlags = NEED_VALUE | NEED_NORMAL | NEED_MEASURE;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, geo));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            geoEval->evaluateAt(quNodes);
            gsMatrix<T> solUVals = loftedVelocity.value(quNodes);
            //gsInfo << "solUVals.col(0) = " << solUVals.col(0) << "\n";

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                gsVector<T> normal;
                geoEval->normal(k, normal);
                normal /= normal.norm();

                const T weight = quWeights(k) * geoEval->measure(k);
                flowRate += weight * solUVals(0, k);// *normal.norm();
            }
        }
        return (m_numBlades * flowRate);
    }

    void computeHead()
    {
        gsWarn << "computeHead() not implemented yet.\n";
    }

    void computeEfficiency()
    {
        gsWarn << "computeEfficiency() not implemented yet.\n";
    }

    void computeCavitation()
    {
        gsWarn << "computeCavitation() not implemented yet.\n";
    }

public:

    void setTargetPressureAtProfile(gsVector<T> pTargetProfile) { m_pTargetProfile = pTargetProfile; }

    void setTargetPressureAtProfile(T shift1, T shift2)
    {
        gsVector<T> pAverage(2);
        for (int k = 0; k < m_profilePatches.rows(); k++)
            pAverage[k] = computeAveragePressure(k);

        T pDiff = math::abs(pAverage[0] - pAverage[1]);
        m_pTargetProfile[0] = pAverage[0] + shift1 * pDiff;
        m_pTargetProfile[1] = pAverage[1] + shift2 * pDiff;
    }

    T computeAveragePressure(int index)
    {
        T averagePressure = 0.;
        T sideLength = 0.;

        const gsBasis<T> & basis = m_bases.at(1).basis(m_profilePatches[index]);

        gsVector<int> numQuadNodes(m_tarDim);
        const int dir = m_profileSides[index].direction();
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = (2 * basis.degree(i) + 1);
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here

        gsMatrix<T> quNodes; // Mapped nodes
        gsVector<T> quWeights; // Mapped weights
        unsigned evFlags = NEED_VALUE | NEED_MEASURE;   // Evaluation flags for the Geometry map

        // Initialize geometry evaluator
          typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_patches.patch(m_profilePatches[index])));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(m_profileSides[index]);
        for (; domIt->good(); domIt->next())
        {
            // Compute the quadrature rule
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            const index_t nQuPoints = quWeights.rows();

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval->evaluateAt(quNodes);

            // Evaluate solution on element nodes
            gsMatrix<T> solPVals = m_pSolField.value(quNodes, m_profilePatches[index]);
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                gsVector<T> normal;
                geoEval->outerNormal(k, m_profileSides[index], normal);

                // weight * abs(det J), where J is geometry Jacobian
                // geoEval->measure(k) and geoEval->jacDet(k) not working as geoEval is initialized at the patch, not on boundary
                const T weight = quWeights[k] * normal.norm();

                averagePressure += weight * solPVals(k);
                sideLength += weight;
            }
        }
        return averagePressure / sideLength;
    }

    /*void setInitialObjectiveParts(T initialLift, T initialOutflowVelocity, T initialPressureProfile)
    {
        m_lift0 = initialLift;
        m_outflowVelocity0 = initialOutflowVelocity;
        m_pressureProfile0 = initialPressureProfile;
    }*/

    void setInitialObjectiveParts(gsVector<T> initialObjectives)
    {
        GISMO_ENSURE(m_tarDim == 2, "Settings of the initial objective parts done only for 2D case.");
        m_lift0 = initialObjectives[0];
        m_outflowVelocity0 = initialObjectives[1];
        m_pressureProfile0 = initialObjectives[2];
        m_objectiveAbsolute = false;
    }

    void setInitialEfficiency(T initialEfficiency)
    {
        m_efficiency0 = initialEfficiency;
        m_objectiveAbsolute = false;
    }

    void setEfficiencyWeight(T weightEfficiency) { m_weightEfficiency = weightEfficiency; }

    void setHead(T head) {m_head = head; }

    void setWeights(T weightLift, T weightOutflowVelocity, T weightPressureProfile, T weightEfficiency = 0.)
    {
        m_weightLift = weightLift;
        m_weightOutflowVelocity = weightOutflowVelocity;
        m_weightPressureProfile = weightPressureProfile;
        m_weightEfficiency = weightEfficiency;

        gsInfo << "m_tarDim = " << m_tarDim << "\n";
        gsInfo << "checkSumOfWeights() = " << checkSumOfWeights() << "\n";

        GISMO_ENSURE(m_tarDim == 2 && checkSumOfWeights() == true, "Dimension or the sum of weights does not match.");
    }

    void setWeightsProfileExample()
    {
        m_weightLift = 0.2;
        m_weightOutflowVelocity = 0.8;
        m_weightPressureProfile = 0.0;
    }

    /*void setWeights(T weightHead, T weightEfficiency, T weightOutflowVelocity, T weightCavitation)
    {
        m_weightHead = weightHead;
        m_weightEfficiency = weightEfficiency;
        m_weightOutflowVelocity = weightOutflowVelocity;
        m_weightCavitation = weightCavitation;

        GISMO_ENSURE(m_tarDim == 3 && checkSumOfWeights() == true, "Dimension or the sum of weights does not match.");
    }*/

    void setWeightsTurbineExample()
    {
        m_weightHead = 0.3;
        m_weightEfficiency = 0.3;
        m_weightOutflowVelocity = 0.3;
        m_weightCavitation = 0.1;
    }

    bool checkSumOfWeights()
    {
        T sumOfWeights = 0.0;
        switch (m_tarDim)
        {
        case 2:
            sumOfWeights = m_weightLift + m_weightOutflowVelocity + m_weightPressureProfile + m_weightEfficiency;
            break;
        case 3:
            sumOfWeights = m_weightHead + m_weightEfficiency + m_weightOutflowVelocity + m_weightCavitation;
            break;
        default:
            GISMO_ERROR("Optimization implemented only for 2D and 3D.");
            break;
        }

        T tol = 1e-15;
        T weightCheck = sumOfWeights - 1.0;
        if (weightCheck < 0)
            weightCheck *= -1;

        if (weightCheck < tol)
            return 1;
        else
            return 0;
    }

public:
    const T getLift() const
    {
        GISMO_ASSERT(m_lift != 0., "Lift equals zero in the denominator.");
        if (isObjectiveAbsolute())
            return m_radius * m_lift;
        else
        {
            if (m_lift0 * m_lift < 0)
                return 1e19;
            else
                return m_lift0 / (m_radius * m_lift);
        }
    }
    const T getInitialLift() const { return m_lift0; }
    const T getOutflowVelocity() const
    {
        if (m_velAvg)
        {
            GISMO_ASSERT(m_uTangentialTarget != 0, "Target velocity equals zero in the denominator.");
            //return m_outflowVelocity / m_uTangentialTarget;
            return m_outflowVelocity;
        }
        else
        {
            GISMO_ASSERT(m_outflowVelocity0 != 0, "Initial velocity objective part equals zero in the denominator.");
            return m_outflowVelocity / m_outflowVelocity0;
        }
    }
    const T getInitialOutflowVelocity() const { return m_outflowVelocity0; }
    const T getTargetVelocity() const { return m_uTangentialTarget; }
    const T getPressureProfile() const
    {
        GISMO_ASSERT(m_pressureProfile0 != 0, "Initial pressure objective part equals zero in the denominator.");
        return m_pressureProfile / m_pressureProfile0;
    }
    const T getInitialPressureProfile() const { return m_pressureProfile0; }
    const gsVector<T> getPressureTargetProfile() const { return m_pTargetProfile; }
    const T getHead() const { return m_head; }
    const T getEfficiency() const { return m_efficiency; }
    const T getCavitation() const { return m_cavitation; }
    const gsVector<T> getObjectives() const
    {
        GISMO_ENSURE(m_objectiveEvaluated == true, "Objective function is not evaluated. Call evaluate() first.");

        gsVector<T> objectives;
        switch (m_tarDim)
        {
        case 2:
            objectives.setZero(3);
            objectives << getLift(), getOutflowVelocity(), getPressureProfile();
            return objectives;
            break;
        case 3:
            objectives.setZero(4);
            objectives << m_head, m_efficiency, m_outflowVelocity, m_cavitation;
            return objectives;
            break;
        default:
            GISMO_ERROR("Optimization implemented only for 2D and 3D.");
            break;
        }
    }

    const bool isObjectiveAbsolute() const { return m_objectiveAbsolute; }
    const index_t getNumObjParts() const
    {
        switch (m_tarDim)
        {
        case 2:
            return 3;
            break;
        case 3:
            return 4;
            break;
        default:
            GISMO_ERROR("Optimization implemented only for 2D and 3D.");
            break;
        }
    }

    const gsVector<T> getWeights() const
    {
        gsVector<T> weights;
        switch (m_tarDim)
        {
        case 2:
            weights.setZero(3);
            weights << m_weightLift, m_weightOutflowVelocity, m_weightPressureProfile;
            return weights;
            break;
        case 3:
            weights.setZero(4);
            weights << m_weightHead, m_weightEfficiency, m_weightOutflowVelocity, m_weightCavitation;
            return weights;
            break;
        default:
            GISMO_ERROR("Optimization implemented only for 2D and 3D.");
            break;
        }
    }

protected:
    T m_weightLift;
    T m_lift, m_lift0;
    gsVector<T> m_normalStream;
    T m_radius;

    T m_weightOutflowVelocity;
    T m_outflowVelocity, m_outflowVelocity0, m_averageVel;
    T m_uTangentialTarget;
    bool m_velAvg;

    T m_weightPressureProfile;
    T m_pressureProfile, m_pressureProfile0;
    gsVector<T> m_pTargetProfile;

    T m_weightHead;
    T m_head;
    T m_targetHead;
    T m_weightEfficiency;
    T m_efficiency, m_efficiency0;
    T m_weightCavitation;
    T m_cavitation;

    int m_tarDim;

    int m_numBlades;
    T m_angularVel;
    T m_flowRate;
    T m_gravity;
    T m_density;

    gsVector<int> m_profilePatches;
    std::vector<boxSide> m_profileSides;
    gsVector<int> m_bladePatches;
    std::vector<boxSide> m_bladeSides;
    gsVector<int> m_outPatches;
    std::vector<boxSide> m_outSides;
    gsVector<int> m_inPatches;
    std::vector<boxSide> m_inSides;

    gsField<T> m_uSolField;
    gsField<T> m_pSolField;
    gsField<T> m_loftedPressure;
    std::vector< gsMultiBasis<T> > m_bases;
    //gsBasis<T>& m_basis;
    gsMultiPatch<T> m_patches;

    bool m_membersReady;
    bool m_objectiveEvaluated;
    bool m_objectiveAbsolute;
};


} // namespace gismo

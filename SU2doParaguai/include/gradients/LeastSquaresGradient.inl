#include "LeastSquaresGradient.hpp"
#include <limits>
#include <cmath>


namespace FVGradient
{
	namespace //just some aliases
	{
		constexpr ulong K11 = 0;
		constexpr ulong K12 = 1;
		constexpr ulong K22 = 2;
	}

	CLeastSquaresGradient::CLeastSquaresGradient(const GradOptions& opts, FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables)
	: CGradient(mesh, variables) 
	{
		ulong nPoin = mesh->GetNumPoints();
		ulong nEdges = mesh->GetNumEdges();
		ulong nDime = mesh->GetDimensions();
		
		invKMat = FVUtil::Block2DMatrix::AllocateBlock2DMatrix(nPoin, 3); 
		//clear invKMat matrix
		for (ulong i = 0; i < nPoin; i++)
		{
			for (uint j = 0; j < 3; j++)
				invKMat(i, j) = 0.0;
		}

		double eps = std::pow(std::numeric_limits<double>::epsilon(), 2); //same eps as su2
		//eps = -1e99;
		weights = new double[nEdges]; // w_ij = w_ji => nEdges weights
		
		//precompute inverse matrix at each node and alternate methods based on determinant
		if (mesh->isMeshStatic()) // this is kinda stupid. I am def not coding a dynamic mesh
		{
			//Acumular deltaX, deltaY e os pesos wij
			//aqui ta tudo armazenado como se fosse invK ja. O determinante é o mesmo, ja que det(invK) * det(K) = 1
			for (ulong e = 0; e < mesh->GetNumEdges(); e++)
			{
				ulong nBeg = mesh->GetEdge(e).nBeg;
				ulong nEnd = mesh->GetEdge(e).nEnd;

				//calcular delta x
				const double delta_x = mesh->GetNode(nEnd).x - mesh->GetNode(nBeg).x;
				//calcular delta y
				const double delta_y = mesh->GetNode(nEnd).y - mesh->GetNode(nBeg).y;
				//calcular peso. Talvez otimizar com outro loop? Teria que guardar os pesos! Dai precisaria pegar delta_x and delta_y de novo. Dane-se isso tambem, nao é simulador da nasa
				const double w_ij_sqr =
					(opts.weighted) 
					? (1.0/(std::pow(delta_x, 2) + std::pow(delta_y, 2))) 
					: 1.0;

				//somar invK11 = w^2 delta_y^2
				const double wdydy = w_ij_sqr * std::pow(delta_y, 2);
				invKMat(nBeg, K11) += wdydy;
				invKMat(nEnd, K11) += wdydy;
				//somar invK12
				const double wdxdy = -w_ij_sqr * delta_y * delta_x; //negativo vem da inversao
				invKMat(nBeg, K12) += wdxdy; 
				invKMat(nEnd, K12) += wdxdy;
				//somar invK22
				const double wdxdx = w_ij_sqr * std::pow(delta_x, 2);
				invKMat(nBeg, K22) += wdxdx;
				invKMat(nEnd, K22) += wdxdx;
				
				weights[e] = w_ij_sqr;
				
			}
			//calcular os determinantes, dividir pelo determinante; se det < eps, invK = 0
			for (ulong i = 0; i < nPoin; i++)
			{
				double &k11 = invKMat(i, K11),
					&k12 = invKMat(i, K12),
					&k22 = invKMat(i, K22);
				//uma hora pra descobrir que eu tinha que ter colocado & atras de cada nome
				const double det = k11 * k22 - k12 * k12;
				if (det < eps)
				{
					k11 = 0; k22 = 0; k12 = 0; //set gradient to zero
				}
				else
				{
					k11 /= det; k12 /= det; k22 /= det; //finalize inverse matrix
				}

			}
			
		}
	}

	void CLeastSquaresGradient::ComputeGradients(const GradOptions& opts, ulong iVarBegin, ulong iVarEnd)
	{
		//basically do the invK * (A^T W^2) b to get the gradient vector
		//iterate through every point and solve the equation there!
		ulong nPoin = Mesh->GetNumPoints();
		ulong nVarPrim = Variables->GetNumPrims();
		ulong nEdges = Mesh->GetNumEdges();

		iVarEnd = (iVarEnd == -1) ? (Variables->GetNumVars() - 1) : iVarEnd;

		//clear gradient
		for (ulong i = 0; i < nPoin; i++)
		{
			for (ulong iVar = iVarBegin; iVar <= iVarEnd; iVar++)
				Variables->SetGradient(i, iVar, { 0.0, 0.0 });
			if (opts.computePrim)
			{
				for (ulong iPrim = 0; iPrim < nVarPrim; iPrim++)
					Variables->SetGradientPrim(i, iPrim, { 0.0, 0.0 });
			}
		}
		//ideia: acumular o vetor c = A^T W^2 b no gradiente e usar para calcular o gradiente sozinho
		//a contribuicao de cada edge é igual pros dois pontos! -delta_x_ji * w_ij^2 * delta_phi_ji = delta_x_ij * w_ji^2 * delta_phi_ij
		for (ulong e = 0; e < nEdges; e++)
		{
			const ulong i = Mesh->GetEdge(e).nBeg, j = Mesh->GetEdge(e).nEnd;

			const double delta_x_ji = Mesh->GetNode(j).x - Mesh->GetNode(i).x;
			const double delta_y_ji = Mesh->GetNode(j).y - Mesh->GetNode(i).y;

			for (ulong iVar = iVarBegin; iVar <= iVarEnd; iVar++)
			{
				const double phi_i = Variables->GetSolution(i, iVar);
				const double phi_j = Variables->GetSolution(j, iVar);
				const double w_delta_phi_ji = weights[e] * (phi_j - phi_i);
				//contributions
				const double delta_c_x = w_delta_phi_ji * delta_x_ji;
				const double delta_c_y = w_delta_phi_ji * delta_y_ji;

				Variables->AddGradient(i, iVar, {delta_c_x, delta_c_y});
				Variables->AddGradient(j, iVar, {delta_c_x, delta_c_y});
			}
			if (opts.computePrim)
			{
				for (ulong iPrim = 0; iPrim < nVarPrim; iPrim++)
				{
					const double phi_i = Variables->GetPrimitive(i, iPrim);
					const double phi_j = Variables->GetPrimitive(j, iPrim);
					const double w_delta_phi_ji = weights[e] * (phi_j - phi_i);
					//contributions
					const double delta_c_x = w_delta_phi_ji * delta_x_ji;
					const double delta_c_y = w_delta_phi_ji * delta_y_ji;
					
					Variables->AddGradientPrim(i, iPrim, {delta_c_x, delta_c_y});
					Variables->AddGradientPrim(j, iPrim, {delta_c_x, delta_c_y});
				}

			}
		}
		
		//usar os valores de invK11,21,22 para achar o gradiente (vetor c dentro do gradiente)
		for (ulong i = 0; i < nPoin; i++)
		{
			for (ulong iVar = iVarBegin; iVar <= iVarEnd; iVar++)
			{
				const double c_x = Variables->GetGradient(i, iVar, FVUtil::COORDINATES::X);
				const double c_y = Variables->GetGradient(i, iVar, FVUtil::COORDINATES::Y);
				
				const double k11 = invKMat(i, K11);
				const double k12 = invKMat(i, K12);
				const double k22 = invKMat(i, K22);

				Variables->SetGradient(i, iVar, 
					{
						c_x * invKMat(i, K11) + c_y * invKMat(i, K12),
						c_x * invKMat(i, K12) + c_y * invKMat(i, K22)
					});
			}
			if (opts.computePrim)
			{
				for (ulong iPrim = 0; iPrim < nVarPrim; iPrim++)
				{
					const double c_x = Variables->GetGradientPrim(i, iPrim, FVUtil::COORDINATES::X);
					const double c_y = Variables->GetGradientPrim(i, iPrim, FVUtil::COORDINATES::Y);

					Variables->SetGradientPrim(i, iPrim, 
						{
							c_x * invKMat(i, K11) + c_y * invKMat(i, K12),
							c_x * invKMat(i, K12) + c_y * invKMat(i, K22)
						});

				}
			}
		}

	}

	

}




#include "TransformCubeToHexagon.h"

void TransformCubeToHexahedron::Transform(std::array<double, SIZE_NODE>& node)
{
	std::array<double, SIZE_NODE> newNode;
	for (int i = 0; i < newNode.size(); i++)
	{
		newNode[i] = 0.0;
		for (int j = 0; j < NUMBER_NODES_CUBE; j++)
			newNode[i] += m_psi(j, node) * m_hexahedron[j][i];
	}
	node = newNode;
}
# FEM-tetrahedrons
В проекте реализован метод конечных элементов на базе элементов-тетраэдров. 

# Сетка
Проект включает в себя построитель сетки, который способен строить сетку, состоящую из элементов-тетраэдров для 
достаточно широкого класса областей. 
## Ограничения и возможности алгоритма построение сетки
1) Расчётная область должны иметь прямолинейные границы. 
2) Внешние границы расчётной области должные составлять выпуклый многогранник. 
3) Расчётная область должна представляться как набор подобластей, которые суть выпуклые шестигранники. 

Также на подобласти и саму расчётную область накладывается ограничение в том, что их нижняя и верхняя границы должны быть параллельны плоскости Oxy.
Границы по осям x и y не обязательно должны быть параллельны осям координат. 

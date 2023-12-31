# Общие сведения
FEM-tetrahedrons - учебный проект, реализующий МКЭ (метод конечных элементов) для решения задач математической физики. 

Используемые технологии:
1. STL
2. RapidJSON
3. Matplotlib (библиотека Python для простой и удобной визуализации результатов)

В проекте реазизован скалярный метод конечных элементов на базе линейных элементов-тетраэдров для эллиптической задачи. 
Основные этапы выполнения и результаты описаны в папке reports.
Изменение параметров программы происходит в res/input. 
Тестирование можно произвести путём указания тестируемой функции в src/main.cpp при создании объекта класса FEM.
Визуализировать конечноэлементную сетку и решение можно в src/visualization.

## Сетка
Проект включает в себя построитель сетки, который аппроксимирует расчётную область элементами-тетраэдрами для достаточно широкого класса областей.

Пример: 

![Image alt](https://github.com/yabaranov/FEM-tetrahedrons/raw/master/res/graph/grid.png)

## Решение эллиптической задачи
Итоговый результат (значения функции в узлах расчётной области) получается в результате решения конечноэлементной разреженной СЛАУ, например, методом BCGSTABLU (src/SLAE), составленной из локальных матриц и векторов конечных элементов (src/FEM).

Некоторые сечения приближения sin(x + y + z) на расчётной области показанной выше (в виде изолиний):

![Image alt](https://github.com/yabaranov/FEM-tetrahedrons/raw/master/res/graph/isolines.png)



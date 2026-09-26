```@meta
CurrentModule = QCLNEGF
```

# [Литература по оптимизации и редукции NEGF](@id optimization-references)

Эта библиография дополняет основную [литературу проекта](references.md) и содержит
только первичные статьи, журнальные описания программных систем и официальные
страницы издателей.  Краткие аннотации ниже описывают **результат авторов в
их задаче**, а не обещание такого же ускорения или точности для reference design.

## Рекурсивные методы и выборочные элементы функций Грина

1. **R. Lake, G. Klimeck, R. C. Bowen, D. Jovanovic**, “Single and multiband
   modeling of quantum electron transport through layered semiconductor
   devices”, *Journal of Applied Physics* **81**, 7845–7869 (1997),
   [doi:10.1063/1.365394](https://doi.org/10.1063/1.365394).
   Локальный слоевой базис приводит задачу Дайсона к последовательности
   связанных блоков; стоимость рекурсии линейна по числу слоёв и кубична по
   размеру поперечного блока.

2. **A. Svizhenko, M. P. Anantram, T. R. Govindan, B. Biegel,
   R. Venugopal**, “Two-dimensional quantum mechanical modeling of
   nanotransistors”, *Journal of Applied Physics* **91**, 2343–2354 (2002),
   [doi:10.1063/1.1432117](https://doi.org/10.1063/1.1432117),
   [авторский текст](https://arxiv.org/abs/cond-mat/0111290).
   Работа выводит рекурсивное вычисление нужных блоков запаздывающей и
   корреляционной функций Грина и связывает его с самосогласованным уравнением
   Пуассона.

3. **L. Lin, C. Yang, J. C. Meza, J. Lu, L. Ying, W. E**, “SelInv—An
   algorithm for selected inversion of a sparse symmetric matrix”,
   *ACM Transactions on Mathematical Software* **37**, 40 (2011),
   [doi:10.1145/1916461.1916464](https://doi.org/10.1145/1916461.1916464).
   Выборочное обращение (selected inversion) вычисляет элементы обратной
   матрицы, соответствующие графу заполнения разреженной
   ``LDL^T``-факторизации, не формируя обратную матрицу целиком.

4. **S. Li, S. Ahmed, G. Klimeck, E. Darve**, “Computing entries of the
   inverse of a sparse matrix using the FIND algorithm”, *Journal of
   Computational Physics* **227**, 9408–9427 (2008),
   [doi:10.1016/j.jcp.2008.06.033](https://doi.org/10.1016/j.jcp.2008.06.033).
   FIND применяет метод вложенных сечений и двунаправленный обход дерева
   разделителей к выбранным элементам функции Грина в широких двумерных и
   трёхмерных устройствах.

5. **D. Mamaluy, D. Vasileska, M. Sabathil, T. Zibold, P. Vogl**,
   “Contact block reduction method for ballistic transport and carrier
   densities of open nanostructures”, *Physical Review B* **71**, 245321 (2005),
   [doi:10.1103/PhysRevB.71.245321](https://doi.org/10.1103/PhysRevB.71.245321).
   Редукция контактных блоков (contact block reduction, CBR) использует
   низкоранговую структуру граничной собственной энергии в открытой
   баллистической задаче. У замкнутой периодической ячейки reference design с полным
   некогерентным рассеянием в SCBA такого контактного блока нет.

Эти методы требуют большого разреженного или блочно-трёхдиагонального
оператора Дайсона. После текущей локальной проекции на ячейку reference design решает
полную матрицу всего лишь
``N_b\times N_b`` для каждой пары ``(E,k)``; поэтому RGF, SelInv и FIND не
являются ускорением существующего ядра. Они становятся кандидатами только
после перехода к координатному представлению, многозонной или широкой
двумерной либо трёхмерной модели.

## FFT, ускорение итераций неподвижной точки и точность арифметики

6. **M. Frigo, S. G. Johnson**, “The Design and Implementation of FFTW3”,
   *Proceedings of the IEEE* **93**, 216–231 (2005),
   [doi:10.1109/JPROC.2004.840301](https://doi.org/10.1109/JPROC.2004.840301).
   FFTW строит аппаратно адаптированные планы дискретного преобразования Фурье
   (DFT). В QCLNEGF быстрое преобразование Фурье (FFT) используется с
   дополнением нулями для линейной, а не циклической тёплицевой свёртки
   дискретного ядра Крамерса—Кронига.

7. **D. G. Anderson**, “Iterative procedures for nonlinear integral
   equations”, *Journal of the ACM* **12**, 547–560 (1965),
   [doi:10.1145/321296.321305](https://doi.org/10.1145/321296.321305).
   Историческая работа вводит многосекантное ускорение итераций неподвижной
   точки. В проекте ускоренная последовательность обязана проверяться
   исходной, а не экстраполированной SCBA-невязкой.

8. **H. F. Walker, P. Ni**, “Anderson Acceleration for Fixed-Point
   Iterations”, *SIAM Journal on Numerical Analysis* **49**, 1715–1735
   (2011),
   [doi:10.1137/10078356X](https://doi.org/10.1137/10078356X).
   Современная формулировка связывает смешивание Андерсона с многосекантной
   квазиньютоновской интерпретацией, подобной GMRES, и обсуждает практическую
   реализацию.

9. **C. G. Broyden**, “A Class of Methods for Solving Nonlinear Simultaneous
   Equations”, *Mathematics of Computation* **19**, 577–593 (1965),
   [doi:10.1090/S0025-5718-1965-0198670-6](https://doi.org/10.1090/S0025-5718-1965-0198670-6).
   Ранговое квазиньютоновское обновление заменяет повторное построение полной
   матрицы Якоби.
   Метод Бройдена не входит в поддерживаемые транспортные алгоритмы.

10. **E. Carson, N. J. Higham**, “Accelerating the Solution of Linear Systems
   by Iterative Refinement in Three Precisions”, *SIAM Journal on Scientific
   Computing* **40**, A817–A847 (2018),
   [doi:10.1137/17M1140819](https://doi.org/10.1137/17M1140819).
   Статья задаёт условия, при которых низкоточная факторизация с
   высокоточной невязкой и итерационным уточнением достигает
   заданной обратной точности.
   Это не разрешение безусловно заменить резонансное решение уравнения Дайсона
   арифметикой FP32.

## Функции Ванье, модальное пространство, LRA и редукция порядка модели

11. **N. Marzari, D. Vanderbilt**, “Maximally localized generalized Wannier
    functions for composite energy bands”, *Physical Review B* **56**,
    12847–12865 (1997),
    [doi:10.1103/PhysRevB.56.12847](https://doi.org/10.1103/PhysRevB.56.12847).
    Унитарная локализация сохраняет выбранное подпространство зон; обрезка
    числа состояний уже вносит погрешность усечения базиса.

12. **S.-C. Lee, A. Wacker**, “Nonequilibrium Green's function theory for
    transport and gain properties of quantum cascade structures”,
    *Physical Review B* **66**, 245314 (2002),
    [doi:10.1103/PhysRevB.66.245314](https://doi.org/10.1103/PhysRevB.66.245314),
    [авторский текст](https://arxiv.org/abs/cond-mat/0212059).
    Метод QCL-NEGF в локализованном базисе включает в SCBA рассеяние на
    фононах, примесях и шероховатости интерфейсов и показывает, почему тензоры
    взаимодействия должны преобразовываться вместе с гамильтонианом.

13. **G. Mil’nikov, N. Mori, Y. Kamakura**, “Equivalent transport models in
    atomistic quantum wires”, *Physical Review B* **85**, 035317 (2012),
    [doi:10.1103/PhysRevB.85.035317](https://doi.org/10.1103/PhysRevB.85.035317).
    Вариационное модальное пространство (variational mode space) строит малый
    представительный базис для заданного энергетического окна и проверяет
    эквивалентность транспортных зон.

14. **D. A. Lemus, J. Charles, T. Kubis**, “Mode-space-compatible inelastic
    scattering in atomistic nonequilibrium Green's function
    implementations”, *Journal of Computational Electronics* **19**,
    1389–1398 (2020),
    [doi:10.1007/s10825-020-01549-8](https://doi.org/10.1007/s10825-020-01549-8),
    [авторский текст](https://arxiv.org/abs/2003.09536).
    Работа рассматривает совместимость неупругого рассеяния с сокращённым
    модальным пространством. Результаты конкретного нанопровода не задают
    автоматически ускорение или допустимое усечение для ядер Фрёлиха ККЛ.

15. **L. Zeng, Y. He, M. Povolotskyi, X. Y. Liu, G. Klimeck, T. Kubis**,
    “Low Rank Approximation Method for Efficient Green's Function Calculation
    of Dissipative Quantum Transport”, *Journal of Applied Physics* **113**,
    213707 (2013),
    [doi:10.1063/1.4809638](https://doi.org/10.1063/1.4809638),
    [авторский текст](https://arxiv.org/abs/1304.0316).
    Глобальное низкоранговое приближение (low-rank approximation, LRA) решает
    NEGF в сокращённом квазичастичном подпространстве и сообщает
    ускорение до ``150\times`` в рассмотренных устройствах. Обратное
    преобразование может давать пространственные осцилляции локального тока,
    поэтому требуются отдельные проверки законов сохранения и погрешности.

16. **J. Z. Huang, W. C. Chew, J. Peng, C.-Y. Yam, L. J. Jiang,
    G. H. Chen**, “Model order reduction for multiband quantum transport
    simulations and its application to p-type junctionless transistors”,
    *IEEE Transactions on Electron Devices* **60**, 2111–2119 (2013),
    [doi:10.1109/TED.2013.2260546](https://doi.org/10.1109/TED.2013.2260546).
    Поперечные моды Блоха проектируют многозонный ``k\cdot p``-гамильтониан в
    существенно меньшее транспортное пространство.

17. **Q. Chen, J. Li, C. Yam, Y. Zhang, N. Wong, G. Chen**, “An approximate
    framework for quantum transport calculation with model order reduction”,
    *Journal of Computational Physics* **286**, 49–61 (2015),
    [doi:10.1016/j.jcp.2015.01.032](https://doi.org/10.1016/j.jcp.2015.01.032),
    [авторский текст](https://arxiv.org/abs/1411.0792).
    Нелинейная редукция порядка модели на основе проекции сокращает число
    дорогих решений по энергии; это суррогатная интерполяция в пространстве
    параметров и энергии, а не тождественная перестановка исходной квадратуры.

Реализованный в QCLNEGF `contraction: low_rank` отличается от работ
11–17: SVD применяется к **готовой развёрнутой матрице дискретного ядра**, а
гамильтониан, функции Грина и число локальных состояний ячейки не
проектируются.

## Приближения микроскопической собственной энергии и другие модели транспорта

18. **T. Kubis, P. Vogl**, “Assessment of approximations in nonequilibrium
    Green's function theory”, *Physical Review B* **83**, 195304 (2011),
    [doi:10.1103/PhysRevB.83.195304](https://doi.org/10.1103/PhysRevB.83.195304).
    Систематически сравниваются локальное, усечённое, усреднённое по импульсу
    и развязанное приближения собственной энергии. Это основной источник для ветвей
    `self_energy_structure` и `transverse_momentum`.

19. **T. Grange**, “Contrasting influence of charged impurities on transport
    and gain in terahertz quantum cascade lasers”, *Physical Review B* **92**,
    241306(R) (2015),
    [doi:10.1103/PhysRevB.92.241306](https://doi.org/10.1103/PhysRevB.92.241306).
    В исследованном терагерцовом QCL ток слабо, а усиление сильно зависит от
    импульсно-зависимого рассеяния на примесях; совпадения только кривой I–V
    недостаточно для оправдания усреднения по импульсу.

20. **M. Büttiker**, “Role of quantum coherence in series resistors”,
    *Physical Review B* **33**, 3020–3026 (1986),
    [doi:10.1103/PhysRevB.33.3020](https://doi.org/10.1103/PhysRevB.33.3020).
    Фиктивный резервуар-зонд заменяет неупругое событие условием нулевого
    полного тока через зонд и позволяет моделировать дефазировку через упругое
    решение.

21. **P. Greck, S. Birner, B. Huber, P. Vogl**, “Efficient method for the
    calculation of dissipative quantum transport in quantum cascade lasers”,
    *Optics Express* **23**, 6587–6600 (2015),
    [doi:10.1364/OE.23.006587](https://doi.org/10.1364/OE.23.006587).
    Метод Бюттикера с множественным рассеянием заменяет квазиравновесным
    выражением меньшую собственную энергию и потому на порядки дешевле полной
    SCBA, но решает другую модель транспорта.

22. **A. Pan, B. A. Burnett, C. O. Chui, B. S. Williams**, “Density matrix
    modeling of quantum cascade lasers without an artificially localized
    basis: A generalized scattering approach”, *Physical Review B* **96**,
    085308 (2017),
    [doi:10.1103/PhysRevB.96.085308](https://doi.org/10.1103/PhysRevB.96.085308).
    Обобщённый супероператор рассеяния переносит микроскопическое рассеяние в
    уравнение для матрицы плотности; это самостоятельная марковская модель, а не
    реализация решения Дайсона—Келдыша.

23. **S. Soleimanikahnoj, O. Jonasson, F. Karimi, I. Knezevic**,
    “Numerically efficient density-matrix technique for modeling electronic
    transport in mid-infrared quantum cascade lasers”, *Journal of
    Computational Electronics* **20**, 280–309 (2021),
    [doi:10.1007/s10825-020-01627-x](https://doi.org/10.1007/s10825-020-01627-x),
    [авторский текст](https://arxiv.org/abs/1710.08870).
    Марковская матрица плотности с сохранением положительности воспроизводит
    выбранные результаты NEGF и эксперимента для средневолнового ИК-QCL до
    порога генерации, однако это не универсальная гарантия для терагерцового
    QCL или ширины линии усиления.

24. **Q. Zhang, M. Tang, L. Wang, K. Xia, Y. Ke**, “Random nonequilibrium
    Green's function method for large-scale quantum transport simulation”,
    *Physical Review B* **110**, 155430 (2024),
    [doi:10.1103/PhysRevB.110.155430](https://doi.org/10.1103/PhysRevB.110.155430).
    ``N_s`` случайных суперпозиций исходных мод заменяют полный базис
    источников; показано стандартное отклонение порядка процента даже при
    малом ``N_s`` благодаря самоусреднению большой поперечной системы.

## Параллельные CPU/GPU реализации

25. **S. Steiger, M. Povolotskyi, H.-H. Park, T. Kubis, G. Klimeck**,
    “NEMO5: A Parallel Multiscale Nanoelectronics Modeling Tool”,
    *IEEE Transactions on Nanotechnology* **10**, 1464–1474 (2011),
    [doi:10.1109/TNANO.2011.2166164](https://doi.org/10.1109/TNANO.2011.2166164).
    Многоуровневое параллельное разбиение применяется к атомистическим,
    многофизическим и NEGF-задачам; эффективность зависит от достаточно
    крупных независимых единиц работы.

26. **S. S. Sawant, F. Léonard, Z. Yao, A. Nonaka**,
    “ELEQTRONeX: A GPU-accelerated exascale framework for non-equilibrium
    quantum transport in nanomaterials”, *npj Computational Materials* **11**,
    110 (2025),
    [doi:10.1038/s41524-025-01604-7](https://doi.org/10.1038/s41524-025-01604-7).
    Среда MPI/GPU демонстрирует самосогласованный расчёт NEGF и электростатики
    на числе до 512 GPU для крупных трёхмерных задач и углеродных нанотрубок.
    В самой статье расширения рассеяния оставлены будущей работой, поэтому её
    масштабирование нельзя считать ориентиром для полной SCBA в QCL.

GPU и распределённое исполнение не меняют физику, пока вычисляется тот же
дискретный оператор и проходят физические проверки в арифметике FP64. Однако
малые ``N_b\times N_b`` блоки Дайсона QCLNEGF сами по себе слишком малы для
эффективной загрузки GPU: потенциальный кандидат — пакетная свёртка ядра над
множеством ``(E,k)``, а не обращение одиночной матрицы.

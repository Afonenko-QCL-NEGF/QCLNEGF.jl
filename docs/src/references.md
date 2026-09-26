# [Основная физическая литература](@id physics-references)

Для последовательного чтения сначала полезен обзор Jirauschek–Kubis:
он связывает квантование, рассеяние и разные уровни транспортной модели.
Работа Lee–Wacker раскрывает конкретное матричное NEGF-представление.
Работа с DOI 10.1063/1.5110305 задаёт экспериментальную структуру, а не исчерпывающую
спецификацию всех приближений данного программного пакета.

1. “Thermoelectrically cooled THz quantum cascade laser
   operating up to 210 K”, *Applied Physics Letters* **115**, 010601 (2019),
   [doi:10.1063/1.5110305](https://doi.org/10.1063/1.5110305),
   [авторский текст](https://arxiv.org/abs/1911.06582).
   Геометрия и листовое легирование: описание лучшей структуры перед рис. 3;
   сопоставление профиля и занятого спектра: рис. 3(a).
2. C. Jirauschek and T. Kubis, “Modeling techniques for quantum cascade
   lasers”, *Applied Physics Reviews* **1**, 011307 (2014),
   [doi:10.1063/1.4863665](https://doi.org/10.1063/1.4863665),
   [arXiv:1412.3563](https://arxiv.org/abs/1412.3563).
   Разделы об электронной структуре, механизмах рассеяния и неравновесных
   функциях Грина дают контекст глав [1–12](@ref theory-model).
3. S.-C. Lee and A. Wacker, “Nonequilibrium Green's function theory for
   transport and gain properties of quantum cascade structures”,
   *Physical Review B* **66**, 245314 (2002),
   [doi:10.1103/PhysRevB.66.245314](https://doi.org/10.1103/PhysRevB.66.245314),
   [arXiv:cond-mat/0212059](https://arxiv.org/abs/cond-mat/0212059).
   Матричные функции Грина, борновское рассеяние и линейный оптический отклик.
4. M. Franckié *et al.*, “Two-well quantum cascade laser optimization by
   non-equilibrium Green's function modelling”, *Applied Physics Letters*
   **112**, 021104 (2018),
   [arXiv:1709.09563](https://arxiv.org/abs/1709.09563).
5. T. Kubis and P. Vogl, “Assessment of approximations in nonequilibrium
   Green's function theory”, *Physical Review B* **83**, 195304 (2011),
   [doi:10.1103/PhysRevB.83.195304](https://doi.org/10.1103/PhysRevB.83.195304).
   Отдельное сравнение допущений собственной энергии и модели контактов.
6. A. Wacker, M. Franckié and D. O. Winge, “Nonequilibrium Green's Function
   Model for Simulation of Quantum Cascade Laser Devices Under Operating
   Conditions”, *IEEE JSTQE* **19**, 1200611 (2013),
   [doi:10.1109/JSTQE.2013.2239613](https://doi.org/10.1109/JSTQE.2013.2239613).
7. A. Wacker, “Gain in quantum cascade lasers and superlattices: A quantum
   transport theory”, *Physical Review B* **66**, 085326 (2002),
   [doi:10.1103/PhysRevB.66.085326](https://doi.org/10.1103/PhysRevB.66.085326).

8. N. Marzari and D. Vanderbilt, “Maximally localized generalized Wannier
   functions for composite energy bands”, *Physical Review B* **56**,
   12847–12865 (1997),
   [doi:10.1103/PhysRevB.56.12847](https://doi.org/10.1103/PhysRevB.56.12847),
   [авторский текст](https://arxiv.org/abs/cond-mat/9707145).
   Унитарная калибровка составного зонного подпространства и пространственный
   разброс; локализация сохраняет подпространство до его усечения.

Данная документация выводит собственные соглашения явно; совпадение названия
метода с литературным не устанавливает эквивалентности всех замыканий.
[Библиография численных методов](@ref optimization-references) дополняет
основные физические источники. Воспроизводимые рисунки используют данные
кода и отдельные подписанные модельные примеры, не копии журнальных рисунков.

Документация программных зависимостей:
[Unitful](https://juliaphysics.github.io/Unitful.jl/stable/),
[Documenter](https://documenter.juliadocs.org/stable/),
[FFTW](https://juliamath.github.io/FFTW.jl/stable/),
[CairoMakie](https://docs.makie.org/stable/documentation/backends/cairomakie/),
[HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/).

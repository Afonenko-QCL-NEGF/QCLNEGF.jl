# QCLNEGF.jl

Библиотека на Julia **1.13.x** для расчёта переноса в слоистых квантово-каскадных
гетероструктурах методом неравновесных функций Грина. Версия **0.2.0** включает
типы для размерных параметров, независимые эталонные и оптимизированные операторы,
решатели SCBA и уравнения Пуассона, расчёт наблюдаемых и оптического отклика.

Для воспроизведения отдельного расчёта используйте его точную версию Julia
и сохранённый `Manifest.toml`.

Библиотека решает явно заданную задачу в памяти. Исследования на основе YAML,
сохранение результатов и контрольных точек, отчёты, графики и командный интерфейс
предоставляет
[QCLNEGFRunner.jl](https://github.com/Afonenko-QCL-NEGF/qcl-negf/tree/main/runner).

## Использование

Установите зависимости из каталога рабочей копии.

```console
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

В небольшом прямом расчёте используется учебный набор параметров.

```julia
using QCLNEGF

problem = build_problem(
    physical = reference_parameters(),
    numerical = tutorial_numerics(),
    scattering = ScatteringOptions(LO=false, acoustic=true,
        impurity=false, IFR=false, alloy=false),
)
solution = solve(problem; options=tutorial_options())
solution.status
```

Перед интерпретацией наблюдаемых проверьте статус решения, невязки и физические
критерии. Сеточную сходимость и согласие с экспериментом устанавливают отдельными
проверками. Малая учебная сетка служит для демонстрации уравнений.

Геометрия эталонной структуры взята из статьи «Thermoelectrically cooled THz
quantum cascade laser operating up to 210 K».
Произвольную последовательность слоёв можно задать через `PhysicalParameters`.
[DOI 10.1063/1.5110305](https://doi.org/10.1063/1.5110305)

## Разработка

```console
deno task check
deno task test
deno task docs
```

Окружение разработки задают проверенный манифест зависимостей и точная версия
Julia. Численные тесты проверяют границы массивов, независимые аналитические
примеры, законы сохранения и согласие прямых и оптимизированных операторов.
Руководство по настройке лаборатории находится в
[qcl-negf/lab](https://github.com/Afonenko-QCL-NEGF/qcl-negf/tree/main/lab).

Документацию для браузера собирает `deno task docs`. Начните с описания
[физической модели](docs/src/theory/01_model.md),
[расчёта в памяти](docs/src/theory/17_implementation_plan.md),
[справочника операторов](docs/src/api/operators.md) и
[архитектуры](docs/src/developer/architecture.md).

Лицензия MIT.

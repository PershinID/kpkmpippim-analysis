Набор скриптов на основе фреймворка ROOT (https://root.cern/) для анализа процесса e+e- аннигиляции e+e- -> K+K-pi+pi-

Папка Base содержит базовые классы:
- для инкапсуляции доступа к деревьям TTree, хранящихся в ROOT-файлах (Container)
- для доступа к одномерным данным и многомерным данными с изменяющимися размерами (Variable), а также расчёта переменных по требованию (TriggerVariable)
- для настройки цикла, в котором происходит расчёт переменных в случае доступа к ним, отбор данных и построение гистограмм (Analysis)

Папка Containers содержит несколько реализаций класса Container с разными наборами переменных, извлекаемых из файлов
Папка Analyses содержит несколько реализаций класса Analysis с разными наборами переменных, которые доступны для анализа или которые требуется рассчитать
Папка Scripts содержит скрипты на базе ROOT и базовых классов

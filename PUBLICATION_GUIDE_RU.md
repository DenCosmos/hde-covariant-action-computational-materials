# Загрузка V52-R1 в GitHub

Инструкция к архиву `V52_R1_GitHub_ready_All_Rights_Reserved.zip`. Научная редакция — V52, вычислительный комплект — V52-R1. Проверка приведённых официальных инструкций: 27 сентября 2026 года.

**Сейчас используется только GitHub. Zenodo, резервирование DOI и его интеграция с GitHub не требуются.** Никаких удалённых действий при подготовке этого архива не выполнено.

## 1. Распаковать подготовленный комплект

Распакуйте архив. Внутри находится папка `V52_R1_repository`. В её корне уже подготовлены `LICENSE`, `LICENSE_SCOPE.md`, `README.md`, `CITATION.cff`, `CITATION.bib`, `SHA256SUMS` и рабочие каталоги.

В `LICENSE` внесён согласованный текст: все права сохраняются за правообладателем в пределах действующего права, с ограниченным разрешением скачивать, хранить, просматривать и запускать неизменённые материалы для личной академической некоммерческой проверки результатов. Имя `Danylo Yerokhin` соответствует метаданным статьи. Права на сторонние библиотеки этим текстом не изменяются.

Перед публичной загрузкой прочитайте [LICENSE_SCOPE.md](LICENSE_SCOPE.md). Собственная лицензия не отменяет условия GitHub. Они предусматривают создание форков средствами сервиса, а также права GitHub и его аффилированных лиц на использование материалов для разработки и обучения ИИ. Отказ от обучения на вводе и выводе отдельных ИИ-функций не отменяет соответствующие положения для публичного репозитория. Решение принять эти условия при размещении материалов остаётся за вами. [G1]

## 2. Создать пустой репозиторий

Войдите в GitHub. Откройте меню **+ → New repository**. Выберите свой аккаунт и имя репозитория. Возможное имя: `hde-covariant-action-materials`; это предложение, а не уже созданный адрес.

В описание можно вставить:

> Computational materials and mathematical supplement for the V52 article on a minimal local covariant action for holographic dark energy. Version V52-R1. All rights reserved; limited scholarly verification permission.

Выберите **Public** только после принятия условий доступа. Не включайте автоматическое добавление README и `.gitignore`. В выборе лицензии оставьте **No license**: это означает не добавлять шаблон GitHub, поскольку ваш собственный файл LICENSE уже подготовлен. Нажмите **Create repository**. [G2, G3]

Результат — пустой репозиторий с настоящим адресом. Сохраните этот адрес; он пригодится для библиографической ссылки. Сам этот шаг ещё не загружает научные файлы.

## 3. Перенести файлы, сохранив каталоги

Для данного вложенного комплекта удобнее GitHub Desktop. В приложении откройте **File → Clone repository**, выберите созданный репозиторий либо вкладку **URL**, задайте локальную папку и нажмите **Clone**. [G4]

Скопируйте в эту локальную папку **содержимое** `V52_R1_repository`, а не ZIP и не саму внешнюю папку. Каталог `.git`, созданный GitHub Desktop, не удаляйте. В корне должны быть непосредственно README, LICENSE и каталоги `scripts`, `lib`, `supplement`, `reference_results`.

Не переносите полный рабочий архив V52, основную статью, окружение Python, пароли и токены. Папка `results/` предназначена для новых локальных запусков и исключена из публикации. Эталонные результаты уже находятся в `reference_results/`.

## 4. Заполнить действительные реквизиты

В `publication/metadata_pending.json` впишите настоящий `repository_url` и выбранную дату предстоящего выпуска `actual_release_date` в формате `YYYY-MM-DD`. После личной проверки авторства и согласия на публичный доступ смените `author_metadata_confirmed` и `public_access_confirmed` на `true`.

Решение о лицензии уже записано: `licenses_confirmed` равно `true`, три области ссылаются на один файл LICENSE. `publication_target` равно `github`. DOI и `license_url` пока могут оставаться `null`. Адрес лицензии добавляется после появления действительной публичной ссылки на этот файл. Не вписывайте вместо него выдуманный адрес.

При использовании Python выполните из корня локального репозитория:

```console
python tools/finalize_metadata.py
python tools/finalize_metadata.py --apply
python tools/validate_citation.py
python tools/check_release.py --update-checksums --for-publication
```

Первый вызов проверяет поля без записи. Второй обновляет метаданные локально, не создавая репозиторий, DOI или релиз. Для маршрута `github` DOI не требуется. Отсутствующие подтверждения или настоящий URL дают `BLOCKED`, а не скрытый пропуск. Команды не запускают научные расчёты. Для работы этих инструментов нужен PyYAML из подготовленного окружения; установку см. в [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

CITATION.cff не содержит фиктивного SPDX-обозначения для нестандартной лицензии. Условия указаны в поле `message` и файле LICENSE. Отсутствие необязательных `license`, `license-url` или `doi` не означает назначения открытой лицензии. Локальный инструмент проверяет YAML, обязательные поля и идентификаторы, а не полную официальную схему CFF. [C1]

## 5. Загрузить содержимое

Во вкладке **Changes** GitHub Desktop проверьте состав. Введите сообщение, например `Prepare V52-R1 materials with restrictive licence`. Нажмите **Commit to main** или кнопку с фактическим именем ветви, затем **Push origin**. Именно Push загружает файлы на сервер. [G5]

В браузере проверьте страницу репозитория, LICENSE, README и PDF дополнения. Скачайте опубликованную копию и сравните `SHA256SUMS` командой `python tools/check_release.py` из корня распакованного репозитория. Название внешней папки скачанного архива может отличаться от `V52_R1_repository`.

Вариант без GitHub Desktop: **Add file → Upload files**. Загружайте распакованные файлы и каталоги партиями с сохранением вложенности. Через браузер действуют ограничения 100 файлов за раз и 25 MiB на файл; данный комплект не помещается в одну партию. Загрузка одного ZIP не заменяет размещение дерева исходников. [G6]

## 6. Зафиксировать версию для статьи

Для ссылки на определённый состав материалов используйте отдельный выпуск. На странице GitHub откройте **Releases → Draft a new release**. Создайте тег `v52-r1`, укажите ветвь с проверенным коммитом и заголовок `V52-R1 computational materials`. В описание внесите научную редакцию V52, условия LICENSE и реальные ограничения из KNOWN_LIMITATIONS.md. DOI не добавляйте, пока его нет. Публикацию выполняете вы кнопкой **Publish release**. [G7]

Можно прикрепить дополнительно ZIP проверенного дерева. Для его сборки, без создания Zenodo-архива, предусмотрена команда:

```console
python tools/build_archives.py --github-only --output-dir ../release-files
```

Архив и `release_record.json` создаются снаружи репозитория. Финальный хеш коммита запишите во внешний протокол и описание релиза, а не внутрь самого идентифицируемого коммита. В статье укажите автора, название материалов, V52-R1 и реальный адрес релиза. Не перемещайте уже процитированный тег на другой коммит.

Zenodo можно подготовить позднее отдельным этапом. Существовавший ранее Zenodo-архив не обновлён при этом лицензионном изменении; его нельзя считать копией текущего GitHub-комплекта.

## Официальные источники

[G1] GitHub Terms of Service, D.3–D.5 и J.3: https://docs.github.com/en/site-policy/github-terms/github-terms-of-service

[G2] Creating a new repository: https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository

[G3] Licensing a repository: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository

[G4] Cloning repositories in GitHub Desktop: https://docs.github.com/en/desktop/adding-and-cloning-repositories/cloning-and-forking-repositories-from-github-desktop

[G5] Committing and reviewing changes: https://docs.github.com/en/desktop/making-changes-in-a-branch/committing-and-reviewing-changes-to-your-project-in-github-desktop

[G6] Adding a file to a repository: https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository

[G7] Managing releases: https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository

[C1] Citation File Format schema guide: https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md

We use Google BigQuery to query the public PyPI download statistics, based on [these instructions](https://packaging.python.org/en/latest/guides/analyzing-pypi-package-downloads/).

### Total number of downloads

```sql
SELECT COUNT(*) AS num_downloads
FROM `bigquery-public-data.pypi.file_downloads`
WHERE file.project = 'text2term'
  AND DATE(timestamp)
    BETWEEN DATE_SUB(CURRENT_DATE(), INTERVAL 36 MONTH)
    AND CURRENT_DATE()
```

### Number of monthly downloads since first release
```sql
SELECT
  COUNT(*) AS num_downloads,
  DATE_TRUNC(DATE(timestamp), MONTH) AS `month`
FROM `bigquery-public-data.pypi.file_downloads`
WHERE
  file.project = 'text2term'
  AND DATE(timestamp)
    BETWEEN DATE_TRUNC(DATE_SUB(CURRENT_DATE(), INTERVAL 36 MONTH), MONTH)
    AND CURRENT_DATE()
GROUP BY `month`
ORDER BY `month` DESC
```

### Number of monthly downloads/installations through **pip** installer since first release
```sql
SELECT
  COUNT(*) AS num_downloads,
  DATE_TRUNC(DATE(timestamp), MONTH) AS `month`
FROM `bigquery-public-data.pypi.file_downloads`
WHERE
  file.project = 'text2term'
  AND details.installer.name = 'pip'
  AND DATE(timestamp)
    BETWEEN DATE_TRUNC(DATE_SUB(CURRENT_DATE(), INTERVAL 36 MONTH), MONTH)
    AND CURRENT_DATE()
GROUP BY `month`
ORDER BY `month` DESC
```
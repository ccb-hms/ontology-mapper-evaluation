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
Result [as of April 17, 2024] — `16,737` downloads

### Number of downloads per month since first release
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

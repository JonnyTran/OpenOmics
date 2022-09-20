import os
import sys
from argparse import ArgumentParser
from pathlib import Path

try:
    from pyspark import SparkConf
    from pyspark.pandas import DataFrame as SparkDataFrame
    from pyspark.sql import SparkSession
    from pyspark.sql import functions as F
    from pyspark.sql.types import *
except ImportError:
    print("Must install `pyspark` with `conda/pip install pyspark` if using read_uniprot_xml(). ")


def read_uniprot_xml(filepath: Path, index_col='accession') \
    -> SparkDataFrame:
    spark = SparkSession.getActiveSession()

    sdf = spark.read.format("com.databricks.spark.xml") \
        .option("rowTag", "entry").option("samplingRatio", 0.01).option("excludeAttribute", True) \
        .option("inferSchema", True).option("nullValue", "") \
        .load(f"file://{filepath}")

    sdf = sdf.withColumn("accession", F.col("accession").getItem(0))
    sdf = sdf.withColumn("organism", F.col("organism").getItem('name'))
    sdf = sdf.withColumn("gene", F.col("gene").getItem('name').getItem(0))
    sdf = sdf.withColumn("geneLocation", F.col("geneLocation").getItem('name').getItem(0))

    # TODO extract `protein` and `feature`
    sdf = sdf.drop('organism').drop('organismHost').drop('comment').drop('dbReference').drop('reference') \
        .drop('evidence').drop('protein').drop('proteinExistence').drop('feature')

    print(sdf.printSchema())
    df: SparkDataFrame = sdf.to_pandas_on_spark(index_col=index_col)

    return df


def start_sparksession():
    os.environ["PYSPARK_PYTHON"] = sys.executable
    os.environ['PYSPARK_SUBMIT_ARGS'] = '--packages com.databricks:spark-xml_2.12:0.15.0 ' \
                                        f'--driver-memory {args["driver-memory"]} ' \
                                        f'--executor-memory {args["executor-memory"]} ' \
                                        'pyspark-shell'
    os.environ['PYARROW_IGNORE_TIMEZONE'] = '1'
    spark = SparkSession.builder \
        .appName('OpenOmics') \
        .getOrCreate()
    # spark.sparkContext.setSystemProperty('spark.driver.maxResultSize', '128g')
    spark.conf.set("spark.sql.execution.arrow.pyspark.enabled", "true")
    return spark


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--filepath', '-f', type=str)
    parser.add_argument('--output', '-o', type=str)
    parser.add_argument('--index', type=str, default="accession")
    parser.add_argument('--driver-memory', type=str, default="128G")
    parser.add_argument('--executor-memory', type=str, default="4G")

    args = parser.parse_args()

    spark = start_sparksession()

    df = read_uniprot_xml(filepath=os.path.expanduser(args.filepath))

    df.to_parquet(args.output, index_col=args.index, )

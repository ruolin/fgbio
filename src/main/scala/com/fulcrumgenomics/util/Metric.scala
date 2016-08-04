/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

package com.fulcrumgenomics.util

import java.io.Writer
import java.nio.file.Path

import dagr.commons.reflect.{ReflectionUtil, ReflectiveBuilder}
import htsjdk.samtools.util.FormatUtil

import scala.reflect.runtime.{universe => ru}
import scala.util.{Failure, Success}

object Metric {
  val Delimiter = "\t"

  /** Get the names of the arguments in the order they were defined for the type [T]. */
  def names[T <: Metric](implicit tt: ru.TypeTag[T]): Seq[String] = {
    val clazz             = ReflectionUtil.typeTagToClass[T]
    val reflectiveBuilder = new ReflectiveBuilder(clazz)
    reflectiveBuilder.argumentLookup.ordered.map(_.name)
  }

  /** Construct a metric of the given type.  The arguments given should be a sequence of tuples: the name of the field
    * and value as a string. */
  def build[T <: Metric](args: Seq[(String, String)])(implicit tt: ru.TypeTag[T]): T = {
    val clazz             = ReflectionUtil.typeTagToClass[T]
    val reflectiveBuilder = new ReflectiveBuilder(clazz)
    // set the arguments
    args.foreach { case (name, value) =>
      reflectiveBuilder.argumentLookup.forField(name) match {
        case Some(arg) =>
          val argumentValue = ReflectionUtil.constructFromString(arg.argumentType, arg.unitType, value) match {
            case Success(v) => v
            case Failure(thr) => throw thr
          }
          arg.value_=(argumentValue)
        case None =>
          throw new IllegalArgumentException(s"The class '${clazz.getSimpleName}' did not have a field with name '$name'.")
      }
    }
    // build it.  NB: if arguments are missing values, then an exception will be thrown here
    reflectiveBuilder.build()
  }

  /** Reads metrics from a set of lines.  The first line should be the header with the field names.  Each subsequent
    * line should be a single metric. */
  def read[T <: Metric](lines: Iterator[String], source: Option[String] = None)(implicit tt: ru.TypeTag[T]): Seq[T] = {
    if (lines.isEmpty) throw new IllegalArgumentException(s"No header found in metrics" + source.map(" in source: " + _).getOrElse(""))
    val names = lines.next().split(Delimiter)
    if (lines.isEmpty) Seq.empty
    else {
      lines.zipWithIndex.map { case (line, lineNumber) =>
        val values = line.split(Delimiter)
        if (names.length != values.length) {
          throw new IllegalArgumentException(s"Line ${lineNumber+1} did not have the same length (${values.length}) as the header (${names.length}")
        }
        build[T](names.zip(values))
      }.toSeq
    }
  }

  /** Reads metrics from the given path.  The first line should be the header with the field names.  Each subsequent
    * line should be a single metric. */
  def read[T <: Metric](path: Path)(implicit tt: ru.TypeTag[T]): Seq[T] = read[T](Io.readLines(path), Some(path.toString))

  /** Writes a metric to the given writer.  The first line will be a header with the field names.  Each subsequent
    * line is a single metric.
    */
  def write[T <: Metric](metrics: Seq[T], writer: Writer)(implicit tt: ru.TypeTag[T]): Unit = {
    writer.write(names[T].mkString(Delimiter))
    writer.write("\n")
    metrics.foreach { metric =>
      writer.write(metric.values.mkString(Delimiter))
      writer.write("\n")
    }
  }

  /** Writes a metric to the given path.  The first line will be a header with the field names.  Each subsequent
    * line is a single metric.
    */
  def write[T <: Metric](metrics: Seq[T], path: Path)(implicit tt: ru.TypeTag[T]): Unit = {
    val writer = Io.toWriter(path)
    write(metrics, writer)
    writer.close()
  }
}

/**
  * Base trait for metrics.
  *
  * All classes extending this class should be a case class.  By convention, all fields should be lower case with
  * words separated by underscores.
  */
trait Metric extends Product {
  private lazy val reflectiveBuilder = new ReflectiveBuilder(this.getClass)
  private val formatter = new FormatUtil

  /** Get the names of the arguments in the order they were defined. */
  def names: Seq[String] = {
    this.reflectiveBuilder.argumentLookup.ordered.map(_.name)
  }

  /** Get the values of the arguments in the order they were defined. */
  def values: Seq[String] = productIterator.map(formatValues).toSeq

  /** Override this method to customize how values are formatted. */
  protected def formatValues(value: Any): String = formatter.format(value)

  override def toString: String = names.zip(values).toMap.toString
}
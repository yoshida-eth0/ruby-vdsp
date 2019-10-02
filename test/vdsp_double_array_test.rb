require "test_helper"

class VdspDoubleArrayTest < Minitest::Test
  def test_create
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.new(3)
    c = a.to_da

    assert_instance_of(Vdsp::DoubleArray, a)
    assert_instance_of(Vdsp::DoubleArray, b)
    assert_instance_of(Vdsp::DoubleArray, c)
  end

  def test_clone
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = a.clone
    c = a.dup

    assert a!=b
    assert a!=c
    assert b!=c

    assert_equal a.to_a, b.to_a
    assert_equal a.to_a, c.to_a
    assert_equal a.length, b.length
    assert_equal a.length, c.length

    a[0] = 4.0
    assert a[0]!=b[0]
    assert a[0]!=c[0]
  end

  def test_length
    a = Vdsp::DoubleArray.new(3)
    assert_equal 3, a.length
    assert_equal 3, a.size
  end

  def test_each
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])

    yield_values = []
    a.each {|v|
      yield_values << v
    }
    assert_equal a.to_a, yield_values
    assert_equal a.to_a, a.each.to_a
  end

  def test_enumerable
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = a.map{|v| v * 2}

    assert_equal [2.0, 4.0, 6.0], b
  end

  def test_to_a
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    assert_equal [1.0, 2.0, 3.0], a.to_a
  end

  def test_aref
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    assert_equal 1.0, a[0]
    assert_equal 2.0, a[1]
    assert_equal 3.0, a[2]
  end

  def test_aset
    a = Vdsp::DoubleArray.new(3)
    a[0] = 1.0
    a[1] = 2.0
    a[2] = 3.0

    assert_equal 1.0, a[0]
    assert_equal 2.0, a[1]
    assert_equal 3.0, a[2]
  end

  def test_operator_add
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])

    assert_equal [2.0, 4.0, 6.0], (a + a).to_a
    assert_equal [3.0, 4.0, 5.0], (a + 2.0).to_a
    assert_equal [3.0, 4.0, 5.0], (2.0 + a).to_a
  end

  def test_operator_sub
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])

    assert_equal [0.0, 0.0, 0.0], (a - a).to_a
    assert_equal [-1.0, 0.0, 1.0], (a - 2.0).to_a
    assert_equal [1.0, 0.0, -1.0], (2.0 - a).to_a
  end

  def test_operator_mul
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])

    assert_equal [1.0, 4.0, 9.0], (a * a).to_a
    assert_equal [2.0, 4.0, 6.0], (a * 2.0).to_a
    assert_equal [2.0, 4.0, 6.0], (2.0 * a).to_a
  end

  def test_operator_div
    a = Vdsp::DoubleArray.create([1.0, 2.0, 4.0])

    assert_equal [1.0, 1.0, 1.0], (a / a).to_a
    assert_equal [0.5, 1.0, 2.0], (a / 2.0).to_a
    assert_equal [2.0, 1.0, 0.5], (2.0 / a).to_a
  end

  def test_vramp
    # TODO
  end

  def test_vgen
    # TODO
  end

  def test_blkman_window
    flag = Vdsp::FULL_WINDOW
    c = Vdsp::DoubleArray.blkman_window(4, flag)

    assert c[0]<c[1]
    assert c[1]<c[2]
    assert c[2]>c[3]
  end

  def test_hamm_window
    flag = Vdsp::FULL_WINDOW
    c = Vdsp::DoubleArray.hamm_window(4, flag)

    assert c[0]<c[1]
    assert c[1]<c[2]
    assert c[2]>c[3]
  end

  def test_hann_window
    flag = Vdsp::FULL_WINDOW
    c = Vdsp::DoubleArray.hann_window(4, flag)

    assert c[0]<c[1]
    assert c[1]<c[2]
    assert c[2]>c[3]
  end

  def test_vclr
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c.vclr

    assert_equal 0.0, c[0]
    assert_equal 0.0, c[1]
    assert_equal 0.0, c[2]
  end

  def test_vfill
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c.vfill(4.0)

    assert_equal 4.0, c[0]
    assert_equal 4.0, c[1]
    assert_equal 4.0, c[2]
  end

  def test_maxv
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.maxv

    assert_equal 3.0, a
  end

  def test_maxmgv
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.maxmgv

    assert_equal 3.0, a
  end

  def test_minv
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.minv

    assert_equal 1.0, a
  end

  def test_minmgv
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.minmgv

    assert_equal 1.0, a
  end

  def test_meanv
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.meanv

    assert_equal 2.0, a
  end

  def test_meamgv
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.meamgv

    assert_equal 2.0, a
  end

  def test_measqv
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.measqv

    assert_in_delta 4.666666666666667, a
  end

  def test_mvessq
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.mvessq

    assert_in_delta 4.666666666666667, a
  end

  def test_rmsqv
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.rmsqv

    assert_in_delta 2.160246899469287, a
  end

  def test_sve
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.sve

    assert_equal 6.0, a
  end

  def test_svemg
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.svemg

    assert_equal 6.0, a
  end

  def test_svesq
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.svesq

    assert_equal 14.0, a
  end

  def test_sve_svesq
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    sum, sum_of_squares = c.sve_svesq

    assert_equal 6.0, sum
    assert_equal 14.0, sum_of_squares
  end

  def test_svs
    c = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a = c.svs

    assert_equal 14.0, a
  end

  def test_biquad
    # TODO
  end
end

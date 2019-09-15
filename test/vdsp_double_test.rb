require "test_helper"

class VdspDoubleTest < Minitest::Test
  def test_vsadd
    # c[i] = a[i] + b
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = 2.5
    c = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vsadd(a, 0, 1, b, c, 0, 1, a.length)

    assert_equal 3.5, c[0]
    assert_equal 4.5, c[1]
    assert_equal 5.5, c[2]
  end

  def test_vadd
    # c[i] = a[i] + b[i]
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vadd(a, 0, 1, b, 0, 1, c, 0, 1, a.length)

    assert_equal 3.5, c[0]
    assert_equal 5.5, c[1]
    assert_equal 7.5, c[2]
  end

  def test_vsub
    # c[i] = a[i] - b[i]
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vsub(b, 0, 1, a, 0, 1, c, 0, 1, a.length)

    assert_equal (-1.5), c[0]
    assert_equal (-1.5), c[1]
    assert_equal (-1.5), c[2]
  end

  def test_vsmul
    # c[i] = a[i] * b
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = 2.5
    c = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vsmul(a, 0, 1, b, c, 0, 1, a.length)

    assert_equal 2.5, c[0]
    assert_equal 5.0, c[1]
    assert_equal 7.5, c[2]
  end

  def test_vmul
    # c[i] = a[i] * b[i]
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vmul(a, 0, 1, b, 0, 1, c, 0, 1, a.length)

    assert_equal 2.5, c[0]
    assert_equal 7.0, c[1]
    assert_equal 13.5, c[2]
  end

  def test_vsdiv
    # c[i] = a[i] / b
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = 2.5
    c = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vsdiv(a, 0, 1, b, c, 0, 1, a.length)

    assert_equal 0.4, c[0]
    assert_equal 0.8, c[1]
    assert_in_delta 1.2, c[2]
  end

  def test_svdiv
    # c[i] = a / b[i]
    a = 2.5
    b = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::svdiv(a, b, 0, 1, c, 0, 1, b.length)

    assert_equal 2.5, c[0]
    assert_equal 1.25, c[1]
    assert_in_delta 0.8333333333333334, c[2]
  end

  def test_vaddsub
    # o0[i] = i1[i] + i0[i]
    # o1[i] = i1[i] - i0[i]
    i0 = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    i1 = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    o0 = Vdsp::DoubleArray.new(3)
    o1 = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vaddsub(i0, 0, 1, i1, 0, 1, o0, 0, 1, o1, 0, 1, i0.length)

    assert_equal 3.5, o0[0]
    assert_equal 5.5, o0[1]
    assert_equal 7.5, o0[2]

    assert_equal 1.5, o1[0]
    assert_equal 1.5, o1[1]
    assert_equal 1.5, o1[2]
  end

  def test_vdiv
    # c[i] = a[i] / b[i]
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vdiv(b, 0, 1, a, 0, 1, c, 0, 1, a.length)

    assert_equal 0.4, c[0]
    assert_in_delta 0.5714285714285714, c[1]
    assert_in_delta 0.6666666666666666, c[2]
  end

  def test_vasm
    # d[i] = (a[i] + b[i]) * c
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = 2.0
    d = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vasm(a, 0, 1, b, 0, 1, c, d, 0, 1, a.length)

    assert_equal 7.0, d[0]
    assert_equal 11.0, d[1]
    assert_equal 15.0, d[2]
  end

  def test_vam
    # d[i] = (a[i] + b[i]) * c[i]
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vam(a, 0, 1, b, 0, 1, c, 0, 1, d, 0, 1, a.length)

    assert_equal 7.0, d[0]
    assert_equal 16.5, d[1]
    assert_equal 30.0, d[2]
  end

  def test_vsbsm
    # d[i] = (a[i] - b[i]) * c
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = 2.0
    d = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vsbsm(a, 0, 1, b, 0, 1, c, d, 0, 1, a.length)

    assert_equal (-3.0), d[0]
    assert_equal (-3.0), d[1]
    assert_equal (-3.0), d[2]
  end

  def test_vsbm
    # d[i] = (a[i] - b[i]) * c[i]
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vsbm(a, 0, 1, b, 0, 1, c, 0, 1, d, 0, 1, a.length)

    assert_equal (-3.0), d[0]
    assert_equal (-4.5), d[1]
    assert_equal (-6.0), d[2]
  end

  def test_vmsa
    # d[i] = (a[i] * b[i]) + c
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = 2.0
    d = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vmsa(a, 0, 1, b, 0, 1, c, d, 0, 1, a.length)

    assert_equal 4.5, d[0]
    assert_equal 9.0, d[1]
    assert_equal 15.5, d[2]
  end

  def test_vma
    # d[i] = (a[i] * b[i]) + c[i]
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vma(a, 0, 1, b, 0, 1, c, 0, 1, d, 0, 1, a.length)

    assert_equal 4.5, d[0]
    assert_equal 10.0, d[1]
    assert_equal 17.5, d[2]
  end

  def test_vmsb
    # d[i] = (a[i] * b[i]) - c[i]
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vmsb(a, 0, 1, b, 0, 1, c, 0, 1, d, 0, 1, a.length)

    assert_equal 0.5, d[0]
    assert_equal 4.0, d[1]
    assert_equal 9.5, d[2]
  end

  def test_vsmsma
    # e[i] = (a[i] * b) + (c[i] * d)
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = 2.5
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = 3.0
    e = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vsmsma(a, 0, 1, b, c, 0, 1, d, e, 0, 1, a.length)

    assert_equal 8.5, e[0]
    assert_equal 14.0, e[1]
    assert_equal 19.5, e[2]
  end

  def test_vaam
    # e[i] = (a[i] + b[i]) * (c[i] + d[i])
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = Vdsp::DoubleArray.create([3.0, 5.0, 7.0])
    e = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vaam(a, 0, 1, b, 0, 1, c, 0, 1, d, 0, 1, e, 0, 1, a.length)

    assert_equal 17.5, e[0]
    assert_equal 44.0, e[1]
    assert_equal 82.5, e[2]
  end

  def test_vmmsb
    # e[i] = (a[i] * b[i]) - (c[i] * d[i])
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = Vdsp::DoubleArray.create([3.0, 5.0, 7.0])
    e = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vmmsb(a, 0, 1, b, 0, 1, c, 0, 1, d, 0, 1, e, 0, 1, a.length)

    assert_equal (-3.5), e[0]
    assert_equal (-8.0), e[1]
    assert_equal (-14.5), e[2]
  end

  def test_vsbsbm
    # e[i] = (a[i] - b[i]) * (c[i] - d[i])
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = Vdsp::DoubleArray.create([3.0, 5.0, 7.0])
    e = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vsbsbm(a, 0, 1, b, 0, 1, c, 0, 1, d, 0, 1, e, 0, 1, a.length)

    assert_equal 1.5, e[0]
    assert_equal 3.0, e[1]
    assert_equal 4.5, e[2]
  end

  def test_vasbm
    # e[i] = (a[i] + b[i]) * (c[i] - d[i])
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.create([2.5, 3.5, 4.5])
    c = Vdsp::DoubleArray.create([2.0, 3.0, 4.0])
    d = Vdsp::DoubleArray.create([3.0, 5.0, 7.0])
    e = Vdsp::DoubleArray.new(3)
    Vdsp::UnsafeDouble::vasbm(a, 0, 1, b, 0, 1, c, 0, 1, d, 0, 1, e, 0, 1, a.length)

    assert_equal (-3.5), e[0]
    assert_equal (-11.0), e[1]
    assert_equal (-22.5), e[2]
  end

  def test_vramp
    # TODO
  end

  def test_vgen
    # TODO
  end

  def test_vramp
    # TODO
  end

  def test_vrampmul2
    # TODO
  end

  def test_vrampmuladd
    # TODO
  end

  def test_vrampmuladd2
    # TODO
  end

  def test_vgenp
    # TODO
  end

  def test_vtabi
    # TODO
  end

  def test_blkman_window
    c = Vdsp::DoubleArray.new(4)
    flag = Vdsp::FULL_WINDOW
    Vdsp::UnsafeDouble.blkman_window(c, c.length, flag)

    assert c[0]<c[1]
    assert c[1]<c[2]
    assert c[2]>c[3]
  end

  def test_hamm_window
    c = Vdsp::DoubleArray.new(4)
    flag = Vdsp::FULL_WINDOW
    Vdsp::UnsafeDouble.hamm_window(c, c.length, flag)

    assert c[0]<c[1]
    assert c[1]<c[2]
    assert c[2]>c[3]
  end

  def test_hann_window
    c = Vdsp::DoubleArray.new(4)
    flag = Vdsp::FULL_WINDOW
    Vdsp::UnsafeDouble.hann_window(c, c.length, flag)

    assert c[0]<c[1]
    assert c[1]<c[2]
    assert c[2]>c[3]
  end
end

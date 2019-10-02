require "vdsp/version"
require "vdsp/vdsp"

module Vdsp
  class Error < StandardError; end
  # Your code goes here...

  module Biquad
    Coefficient = Struct.new("Coefficient", :b0, :b1, :b2, :a1, :a2)
  end
end

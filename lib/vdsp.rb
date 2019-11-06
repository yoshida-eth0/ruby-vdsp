require "vdsp/version"
require "vdsp/vdsp"

module Vdsp
  class Error < StandardError; end
  # Your code goes here...

  module Array
    def blkman_window(flag=FULL_WINDOW)
      window = self.class.blkman_window(self.length, flag)
      self * window
    end

    def hamm_window(flag=FULL_WINDOW)
      window = self.class.hamm_window(self.length, flag)
      self * window
    end

    def hann_window(flag=FULL_WINDOW)
      window = self.class.hann_window(self.length, flag)
      self * window
    end
  end

  module Biquad
    Coefficient = Struct.new("Coefficient", :b0, :b1, :b2, :a1, :a2)
  end
end

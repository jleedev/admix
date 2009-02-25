class AlleleFrequency < ActiveRecord::Base
  belongs_to :marker

  def freq_comp
    1 - freq
  end
end

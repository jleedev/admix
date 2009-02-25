class Marker < ActiveRecord::Base
  belongs_to :locus_file
  has_many :allele_frequency

  def export
    "#{name}\n0 #{(allele_frequency.map &:freq).join ' '}\n1 #{(allele_frequency.map &:freq_comp).join ' '}"
  end
end

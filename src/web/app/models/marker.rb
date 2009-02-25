class Marker < ActiveRecord::Base
  belongs_to :locus_file
  has_many :allele_frequency

  def export
    "#{name}\n#{(allele_frequency.map &:freq).join ' '}"
  end
end

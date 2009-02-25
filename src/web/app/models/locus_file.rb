class LocusFile < ActiveRecord::Base
  belongs_to :user
  has_many :marker
  has_many :job
  validates_presence_of :name,:population
  validates_numericality_of :population, :greater_than_or_equal_to => 2

  def export
    ((marker.map &:export).join "\n") + "\n"
  end
end

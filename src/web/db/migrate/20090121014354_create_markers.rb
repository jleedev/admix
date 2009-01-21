class CreateMarkers < ActiveRecord::Migration
  def self.up
    create_table :markers do |t|
      t.string(16) :name

      t.timestamps
    end
  end

  def self.down
    drop_table :markers
  end
end
